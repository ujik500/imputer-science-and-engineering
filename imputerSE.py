import argparse

def main():
    parser = argparse.ArgumentParser(description="Description...")
    parser.add_argument("--labels", type=str, default="snp_data/pop_labels.tsv", help="Path to file with population labels for each sample")
    parser.add_argument("--train", type=str, default="snp_data/chr1_train.vcf.gz", help="Path to training data with fully intact SNP data")
    args = parser.parse_args()
    labels_path = args.labels
    train_path = args.train

    train_file = open(train_path, "r")
    lines = train_file.readlines()
    for i in range(20):
        print(lines[i])



if __name__ == "__main__":
    main()