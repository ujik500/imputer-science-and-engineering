import gzip

def truncate_vcf(input_path, output_path, max_snps=25000):
    # Open necessary files
    with gzip.open(input_path, 'rt') as f_in, open(output_path, 'w') as f_out:
        snp_count = 0
        
        for line in f_in:
            # Skip metadata lines
            if line.startswith('##'):
                continue
            
            # Write the header line
            if line.startswith('#'):
                f_out.write(line)
                continue
            
            # Write SNP data lines
            f_out.write(line)
            snp_count += 1
            
            # Stop once we hit the target number of SNPs
            if snp_count >= max_snps:
                break
                
    print(f"Extraction complete! Saved header and {snp_count} SNPs to {output_path}")

if __name__ == "__main__":
    input_file = "snp_data/chr1_train.vcf.gz"
    output_file = "snp_data/chr1_train_medium_truncated.vcf"
    
    truncate_vcf(input_file, output_file, 5000)