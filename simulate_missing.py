import argparse
import random
import gzip

def simulate_missing(percent_missing, num_samples, train_path):
    out = open("missing" + str(percent_missing) + ".vcf", "w")
    
    header_seen = False
    snp_count = 0
    
    with gzip.open(train_path, 'rt') as f:
        for line in f:
            if line.startswith("##"):
                continue
            elif line.startswith("#") and header_seen == False:
                header = line.strip().split("\t")
                length = len(header)
                header_seen = True
            else:
                snp_line = line.strip().split("\t")
                labels = snp_line[:9]
                test_set = snp_line[len(snp_line)-num_samples:]
                for i in range(len(test_set)):
                    # mask out SNPs with the specified probability
                    if random.random() * 100 < percent_missing:
                        test_set[i] = "???"
                final_line = "\t".join(labels) + "\t" + "\t".join(test_set)
                out.write(final_line + "\n")

                snp_count += 1
                if snp_count % 100000 == 0:
                    print(snp_count)
                if snp_count > 25000:
                    print(line.strip().split("\t"))
                    break
    return


simulate_missing(10, 50, "snp_data/chr1_train.vcf.gz")