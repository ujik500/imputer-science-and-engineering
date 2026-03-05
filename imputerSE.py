import argparse
import gzip

def process_ancDNA(ancDNA_path):
    snp_list = []
    try:
        f = open(ancDNA_path)
        lines = f.readlines()
        for line in lines:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            if line[1] == "1":
                snp_list.append(line[0])
            
    except FileNotFoundError:
        print(f"Error: The file {ancDNA_path} was not found.")
    
    return snp_list
    
def process_training(train_path, test_snps):
    header_seen = False
    snp_count = 0
    overlap = 0

    try:
        with gzip.open(train_path, 'rt') as f:
            for line in f:
                # Process each line here
                if line.startswith("##"):
                    continue
                elif line.startswith("#") and header_seen == False:
                    header = line.strip().split("\t")
                    header_seen = True
                else:
                    snp_line = line.strip().split("\t")
                    if snp_line[2] in test_snps: #check to see if 
                        overlap += 1
                    snp_count += 1
                    if snp_count % 100000 == 0:
                        print(snp_count)
                    if snp_count > 25000:
                        print(line.strip().split("\t"))
                        break
                        
                    position = int(snp_line[1])
                    if abs(position - 1000000) < 1000:
                        print(position)
        print(snp_count)
        print(overlap)
    except FileNotFoundError:
        print(f"Error: The file {train_path} was not found.")

def main():
    parser = argparse.ArgumentParser(description="Description...")
    parser.add_argument("--labels", type=str, default="snp_data/pop_labels.tsv", help="Path to file with population labels for each training sample")
    parser.add_argument("--train", type=str, default="snp_data/chr1_train.vcf.gz", help="Path to training data with fully intact SNP data")
    parser.add_argument("--ancDNA", type=str, default="snp_data/anonAncestryDNA.txt", help="Path to a single AncestryDNA raw datafile to impute for")
    args = parser.parse_args()
    labels_path = args.labels
    train_path = args.train
    ancDNA_path = args.ancDNA
    
    test_snps = process_ancDNA(ancDNA_path)
    process_training(train_path, test_snps)
    



if __name__ == "__main__":
    main()