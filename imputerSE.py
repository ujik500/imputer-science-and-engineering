import gzip
import argparse
from tqdm import tqdm
import matplotlib.pyplot as plt

CORE_SIZE = 1000    # num of SNPs to look at during each segment
CONTEXT = 50        # num of SNPs to look at as context hanging off each side
WINDOW = CORE_SIZE + 2 * CONTEXT # total window size
STEP = CORE_SIZE    # offset for where next window begins
MAX_SNPS = 25000    # num of SNPs to impute for
NUM_TEST_SAMPLES = 50   # num of samples to impute for, and exclude from training set

def read_test_file(pred_path):
    with open(pred_path) as f:
        lines = [l.strip().split("\t") for l in f.readlines()]
    return lines

'''
Iterator to step through each segment of the test data. LLMs helped me the most
when writing this function and the following one, due to complex alignment
logic.
'''
def test_segment_iterator(lines):
    for i in range(0, len(lines), CORE_SIZE):
        start = max(0, i - CONTEXT)
        end = i + CORE_SIZE + CONTEXT
        yield lines[start:end]

'''
Iterator to step through each segment of the training data. Separate logic from
test_segment_iterator() since the headers and metadata of the files are different.
'''
def train_segment_iterator(train_path):
    # handle zipped or unzipped files correctly
    if not train_path.endswith(".gz"):
        with open(train_path, 'rb') as f_in, gzip.open(train_path + ".gz", 'wb') as f_out:
            f_out.write(f_in.read())
        train_path += ".gz"
        
    snp_buffer = []
    header_seen = False
    snp_count = 0

    with gzip.open(train_path, "rt") as f:
        for line in f:
            if line.startswith("##"): continue
            if line.startswith("#") and not header_seen:
                header = line.strip().split("\t")
                gt_start = header.index("FORMAT") + 1
                header_seen = True
                continue

            if snp_count >= MAX_SNPS: break

            # Exclude the last NUM_TEST_SAMPLES columns from training, prevent leakage
            snp_buffer.append(line.strip().split("\t")[:len(header)-NUM_TEST_SAMPLES])
            snp_count += 1

            if len(snp_buffer) >= WINDOW:
                yield snp_buffer[:WINDOW], gt_start
                snp_buffer = snp_buffer[STEP:]

        # Yield the leftovers that didn't make a full segment
        if len(snp_buffer) > 0:
            yield snp_buffer, gt_start

'''
This function takes a streamed segment of a vcf file and the column index 
where genotypes begin, and returns a list of haplotype strings.
'''
def segment_to_haplotypes(seg, gt_start):
    num_samples = len(seg[0]) - gt_start
    num_haps = num_samples * 2

    hap_lists = [[] for hap in range(num_haps)]

    for snp in seg:
        genotypes = snp[gt_start:]

        for i, gt in enumerate(genotypes):
            h1 = 2 * i
            h2 = 2 * i + 1

            if gt == "???":
                a, b = "?", "?"
            else:
                try:
                    a, b = gt.split("|")
                except:
                    a, b = "?", "?"

            hap_lists[h1].append(a)
            hap_lists[h2].append(b)

    return ["".join(h) for h in hap_lists]

'''
This function does a greedy search of train haplotypes to find the one with the lowest
Hamming Distance from each test haplotype. Then, it predicts missing test genotypes 
based on the values in its assigned neighbor. It returns a list of imputed test haplotype strings.
'''
def nearest_haplotype_impute(test_haps, train_haps):
    imputed = []
    for th in test_haps:
        best_match = None
        best_dist = 10**9

        for rh in train_haps:
            dist = 0
            # zip automatically stops at the shortest length
            for a, b in zip(th, rh):
                if a == "?": continue
                if a != b: dist += 1

            if dist < best_dist:
                best_dist = dist
                best_match = rh

        new_hap = list(th)
        # only iterate up to the length of the best_match found
        for i in range(len(best_match)):
            if new_hap[i] == "?":
                new_hap[i] = best_match[i]

        imputed.append("".join(new_hap))
    return imputed


def write_predictions(pred_file, segment_index, haplotypes):
    with open(pred_file, "a") as f:
        for h in haplotypes:
            f.write(f"{segment_index}\t{h}\n")

'''
This function computes accuracy values, separated at SNPs with common and rare (MAF < 0.01)
variants. Recall for alternate alleles is also reported in the same way. A histogram of rare and common
recalls across
'''
def test_accuracy(train_path, test_path, pred_path):
    if not train_path.endswith(".gz"):
        train_path += ".gz"
    correct = {}
    totals = {}
    
    rare_correct = {}
    rare_totals = {}
    
    recall_correct = 0
    total = 0
    
    rare_recall_correct = 0
    rare_total = 0
    
    # Stitch segments together by CORE_SIZE to align indexes
    full_pred_haps = [[] for _ in range(NUM_TEST_SAMPLES * 2)]
    with open(pred_path, "r") as pred_file:
        # Use a counter to track which haplotype in the sample set we are on
        for line_num, line in enumerate(pred_file):
            seg_idx, hap_str = line.strip().split("\t")
            seg_idx = int(seg_idx)
            
            # This identifies which of the 100 haplotypes (50 samples * 2) this line belongs to
            h_idx = line_num % (NUM_TEST_SAMPLES * 2)
            
            start_bit = CONTEXT if seg_idx > 0 else 0
            full_pred_haps[h_idx].append(hap_str[start_bit:start_bit+CORE_SIZE])
    
    # Flatten the stitched segments into one string per haplotype
    pred_haplotypes = ["".join(chunks) for chunks in full_pred_haps]
            
    missing = []
    with open(test_path, "r") as test_file:
        for line in test_file:
            missing_here = []
            line = line.strip().split("\t")
            for i in range(len(line)):
                if line[i] == "???":
                    missing_here.append(i)
            missing.append(missing_here)
            
    snp_idx = 0
    header_seen = False
    with gzip.open(train_path, "rt") as train_file:
        for line in train_file:
            if line.startswith("##"): 
                continue
            if line.startswith("#") and not header_seen:
                header = line.strip().split("\t")
                header = header[len(header)-NUM_TEST_SAMPLES:]
                header_seen = True
                continue
            line = line.strip().split("\t")
            
            # filter out non bi-allelic SNPs for accuracy calculations
            try:
                maf = float(line[7].split(";")[1][3:])
            except ValueError:
                snp_idx += 1
                continue
                
            line = line[len(line)-NUM_TEST_SAMPLES:]
            
            # Stop if we've reached the end of what was actually imputed
            if snp_idx >= len(pred_haplotypes[0]): break

            for i in range(len(line)):
                if i in missing[snp_idx]:                                             
                    name = header[i]
                    allele_1 = line[i][0]
                    allele_2 = line[i][2]

                    pred_allele_1 = pred_haplotypes[i * 2][snp_idx]
                    pred_allele_2 = pred_haplotypes[i * 2 + 1][snp_idx]

                    if maf > 0.01:
                        if allele_1 == pred_allele_1:
                            correct[name] = correct.get(name, 0) + 1
                            if pred_allele_1 != "0":
                                recall_correct += 1

                        if allele_2 == pred_allele_2:
                            correct[name] = correct.get(name, 0) + 1
                            if pred_allele_2 != "0":
                                recall_correct += 1

                        totals[name] = totals.get(name, 0) + 2
                        
                        if allele_1 != "0":
                            total += 1
                        if allele_2 != "0":
                            total += 1    
                    else:
                        if allele_1 == pred_allele_1:
                            rare_correct[name] = rare_correct.get(name, 0) + 1
                            if pred_allele_1 != "0":
                                rare_recall_correct += 1

                        if allele_2 == pred_allele_2:
                            rare_correct[name] = rare_correct.get(name, 0) + 1
                            if pred_allele_2 != "0":
                                rare_recall_correct += 1

                        rare_totals[name] = rare_totals.get(name, 0) + 2
                        
                        if allele_1 != "0":
                            rare_total += 1
                        if allele_2 != "0":
                            rare_total += 1 
                        
            
            snp_idx += 1
    
    print(f"Accuracy across {snp_idx} SNPs:")
    for key in correct:
        correct[key] = round(correct[key] / totals[key], 4)
    for key in rare_correct:
        rare_correct[key] = round(rare_correct[key] / rare_totals[key], 4)
        
    
    print("Accuracy for common variants:", sum(correct.values()) / len(correct))
    print("Accuracy for rare variants:", sum(rare_correct.values()) / len(rare_correct))
    
    print("Recall for common variants:", recall_correct / total)
    print("Recall for rare variants:", rare_recall_correct / rare_total)

    # Extract accuracy values from dictionaries
    common_acc_values = list(correct.values())
    rare_acc_values = list(rare_correct.values())

    # Create histogram of rare and common accuracies for each sample
    plt.figure(figsize=(10, 6))
    plt.hist(common_acc_values, bins=20, alpha=0.5, label='Common Variants', 
             color='blue', edgecolor='black')
    plt.hist(rare_acc_values, bins=20, alpha=0.5, label='Rare Variants', 
             color='orange', edgecolor='black')

    # Add labels and styling
    plt.title('Distribution of Prediction Accuracies across Samples')
    plt.xlabel('Accuracy')
    plt.ylabel('Frequency (Number of Samples)')
    plt.legend(loc='upper left')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Ensure labels are not truncated
    plt.tight_layout()

    # Save or show the plot
    plt.savefig('accuracy_distribution.png')

    
'''
Runs the full SNP genotype imputation pipeline.
'''
def run_pipeline(train_path, test_path, pred_path):
    test_lines = read_test_file(test_path)

    test_iter = test_segment_iterator(test_lines)
    train_iter = train_segment_iterator(train_path)

    open(pred_path, "w").close()

    # Capture all training segments and match lengths for the final segment
    all_train = list(train_iter)
    
    for i, ((train_seg, gt_start), test_seg) in tqdm(enumerate(zip(all_train, test_iter)), desc="Imputing on test segments...", total=len(all_train)):
        train_haps = segment_to_haplotypes(train_seg, gt_start)
        test_haps = segment_to_haplotypes(test_seg, gt_start)

        # Align window sizes if needed (essential for the tail end)
        L = min(len(train_haps[0]), len(test_haps[0]))
        train_haps = [h[:L] for h in train_haps]
        test_haps = [h[:L] for h in test_haps]

        imputed = nearest_haplotype_impute(test_haps, train_haps)

        write_predictions(pred_path, i, imputed)

    test_accuracy(train_path, test_path, pred_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ABCDEF")
    parser.add_argument("--train", type=str, default="snp_data/chr1_train_medium_truncated.vcf", help="Path to phased training data with fully intact SNP data")
    parser.add_argument("--test", type=str, default="masked_data/missing10.vcf", help="Path to phased incomplete data to impute for")
    parser.add_argument("--out", type=str, default="predictions.txt", help="Path to output file")
    parser.add_argument("--window_size", type=int, default=1000, help="Size of each sliding window segment, not counting overhangs")
    args = parser.parse_args()
    

    train_path = args.train
    test_path = args.test
    pred_path = args.out
    CORE_SIZE = args.window_size
    WINDOW = CORE_SIZE + 2 * CONTEXT # total window size
    STEP = CORE_SIZE

    run_pipeline(train_path, test_path, pred_path)