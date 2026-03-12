import gzip
import argparse

CORE_SIZE = 1000    # number of SNPs to look at during each segment
CONTEXT = 50        # number of SNPs to look at as context hanging off each side
WINDOW = CORE_SIZE + 2 * CONTEXT # total window size
STEP = CORE_SIZE
MAX_SNPS = 25000    # number of SNPs to impute for
NUM_TEST_SAMPLES = 50   # number of samples to impute for, and exclude from training set

def read_test_file(pred_path):
    with open(pred_path) as f:
        lines = [l.strip().split("\t") for l in f.readlines()]
    return lines


def test_segment_iterator(lines):
    num_segments = len(lines) // CORE_SIZE + 1

    for seg in range(num_segments):

        start = max(0, seg * CORE_SIZE - CONTEXT)
        end = seg * CORE_SIZE + CORE_SIZE + CONTEXT

        yield lines[start:end]


def train_segment_iterator(train_path):
    if not train_path.endswith(".gz"):
        with open(train_path, 'rb') as f_in, gzip.open(train_path + ".gz", 'wb') as f_out:
            f_out.write(f_in.read())
        train_path += ".gz"
    snp_buffer = []
    header_seen = False
    snp_count = 0

    with gzip.open(train_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#") and not header_seen:
                header = line.strip().split("\t")
                gt_start = header.index("FORMAT") + 1
                header_seen = True
                continue

            if snp_count >= MAX_SNPS:
                break

            line = line.strip().split("\t")
            line = line[:len(line) - NUM_TEST_SAMPLES] # remove test samples from training
            snp_buffer.append(line)
            snp_count += 1

            if len(snp_buffer) < WINDOW:
                continue

            yield snp_buffer[:WINDOW], gt_start

            snp_buffer = snp_buffer[STEP:]


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


def nearest_haplotype_impute(test_haps, train_haps):
    imputed = []

    for th in test_haps:
        best_match = None
        best_dist = 10**9

        for rh in train_haps:
            dist = 0

            for a, b in zip(th, rh):
                if a == "?":
                    continue

                if a != b:
                    dist += 1

            if dist < best_dist:
                best_dist = dist
                best_match = rh

        new_hap = list(th)

        for i, allele in enumerate(new_hap):
            if allele == "?":
                new_hap[i] = best_match[i]

        imputed.append("".join(new_hap))

    return imputed


def write_predictions(pred_file, segment_index, haplotypes):
    with open(pred_file, "a") as f:
        for h in haplotypes:
            f.write(f"{segment_index}\t{h}\n")

def test_accuracy():
    pass

def run_pipeline(train_path, test_path, pred_path):
    test_lines = read_test_file(test_path)

    test_iter = test_segment_iterator(test_lines)
    train_iter = train_segment_iterator(train_path)

    open(pred_path, "w").close()

    for i, ((train_seg, gt_start), test_seg) in enumerate(zip(train_iter, test_iter)):
        train_haps = segment_to_haplotypes(train_seg, gt_start)

        test_haps = segment_to_haplotypes(test_seg, gt_start)

        # Align window sizes if needed
        if len(train_haps[0]) > len(test_haps[0]):
            L = len(test_haps[0])
            train_haps = [h[:L] for h in train_haps]

        imputed = nearest_haplotype_impute(test_haps, train_haps)

        write_predictions(pred_path, i, imputed)

        print(f"Segment {i} imputed")


if __name__ == "__main__":

    train_path = "snp_data/chr1_train_super_truncated.vcf"
    test_path = "missing10.vcf"
    pred_path = "predictions.txt"

    run_pipeline(train_path, test_path, pred_path)

'''
def main():
    parser = argparse.ArgumentParser(description="Description...")
    parser.add_argument("--labels", type=str, default="snp_data/pop_labels.tsv", help="Path to file with population labels for each training sample")
    parser.add_argument("--train", type=str, default="snp_data/chr1_train.vcf.gz", help="Path to phased training data with fully intact SNP data")
    parser.add_argument("--pred", type=str, default="missing10.vcf", help="Path to phased incomplete data to impute for")
    args = parser.parse_args()
    labels_path = args.labels
    train_path = args.train
    pred_path = args.pred
    
    segmented_haplotypes = process_pred(pred_path)
    print(segmented_haplotypes[0][0])
    predict_missing(train_path, segmented_haplotypes)

if __name__ == "__main__":
    main()
'''