# Imputer Science and Engineering
My CSE 284 project, which is a tool to perform imputation on incomplete SNP genotype data.

## Project Overview
I am using phased data from chromosome 1 of the 1000Genomes dataset (minus the last 50 samples) as training data to inform my predictions, and the final 50 samples as test data. To simulate incomplete data, I am randomly masking 10% of the SNP genotypes in the final 50 samples, and saving the results in a file called ```missing10.vcf```, which you can't find in the repo, since it is too large. Right now, it only uses the first 25k SNPs, since runtime is still a serious concern during testing. 

## Example Usage
```
python imputerSE.py --train snp_data/chr1_train_medium_truncated.vcf --test masked_data/missing10.vcf --out my_predictions.txt --window_size 1000
```
Both train and test data must be phased, and all of the options will have defaults if none are provided.

## File Explanations
### Python Scripts:
* ```imputerSE.py```: The main imputation pipeline
* ```truncate_vcf.py```: Code used to extract just the first 25k SNPs of the chr1 1000Genomes vcf. This file can also be used to generate even smaller datasets, simply change the ```max_snps``` parameter.
* ```simulate_missing.py```: Code used to simulate missing data from 1000Genomes at the rate of your choice.

### Data:
* ```missing10.vcf```: chr1 1000Genomes data with 10% of SNP genotypes randomly masked out
* ```snp_data/chr1_train_super_truncated.vcf```: The first 1000 SNPs on chr1, taken from 1000Genomes
* ```snp_data/chr1_train_truncated.vcf```: The first 25000 SNPs on chr1, taken from 1000Genomes




Also, one of my goals for the project was being able to predict missing SNP genotypes on my own AncestryDNA dataset. However, the SNP data from Ancestry or 23andMe is very sparse compared to the 1000Genomes training set, so there can be large gaps (multiple kb), so I would guess that the power to impute with LD is greatly reduced. I will have to decide how to approach that, since I also don't have a very accurate way to benchmark it.


