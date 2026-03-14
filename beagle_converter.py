import gzip

f = open("snp_data/chr1_train_super_truncated.vcf", "r")
header = f.readline()
header = header.strip().split("\t")
header = header[:9] + header[len(header)-50:]
header = "\t".join(header)
f.close()

f = open("missing10.vcf", "r")
lines = f.readlines()
f.close()
f = open("beagle_inputs/missing10_beagle.vcf", "w")
f.write(header + "\n")
for line in lines:
    line = line.strip().split("\t")
    if len(line) != 59:
        print(len(line))
    for i in range(len(line)):
        if line[i] == "???" or line[i] == "./.":
            line[i] = ".|."
    line = "\t".join(line)
    f.write(line + "\n")
f.close()

with gzip.open("snp_data/chr1_train.vcf.gz", "rt") as f, open("beagle_inputs/chr1_train_truncated_beagle.vcf", "w") as out:
    counter = 0
    for line in f:
        if not line.startswith("#"):
            counter += 1
            
        if line.startswith("##"):
            continue
        else:
            line = line.strip().split("\t") #.decode("utf-8")
            line = line[:len(line)-50]
            out.write("\t".join(line) + "\n")
        
        if counter == 25000:
            break
f.close()
out.close()