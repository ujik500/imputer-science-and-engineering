# Imputer Science and Engineering
My CSE 284 project, which is a tool to perform imputation on incomplete SNP genotype data.

## Project Overview
I am using phased data from chromosome 1 of the 1000Genomes dataset (minus the last 100 samples) as training data to inform my predictions, and the final 100 samples as test data. To simulate incomplete data, I am randomly masking 10% of the SNP genotypes in the final 100 samples, and saving the results in a file called ```missing10.vcf```, which you can find in the repo. So far, it only uses the first 25k SNPs, since runtime is still a serious concern during testing. Using the command ```python imputerSE.py``` (simple for now, but will convert into a full command-line tool for the final), you can run the same analysis pipeline, which takes the truncated chromosome 1 data, and imputes.

## File Explanations
### Python Scripts:
* ```imputerSE.py```: The main imputation pipeline
* ```truncate_vcf.py```: Code used to extract just the first 25k SNPs of the chr1 1000Genomes vcf. This file can also be used to generate even smaller datasets, simply change the ```max_snps``` parameter.
* ```simulate_missing.py```: Code used to simulate missing data from 1000Genomes at the rate of your choice.

### Data:
* ```missing10.vcf```: chr1 1000Genomes data with 10% of SNP genotypes randomly masked out
* ```snp_data/chr1_train_super_truncated.vcf```: The first 1000 SNPs on chr1, taken from 1000Genomes
* ```snp_data/chr1_train_truncated.vcf```: The first 25000 SNPs on chr1, taken from 1000Genomes

Other files can be ignored for now.

## Other Notes for Peer Reviewers
I would recommend that you clone this repo in your CSE 284 Datahub, so you can easily access the data and compute resources needed. You will need to move ```~/public/1000Genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz``` into this repo's folder ```snp_data``` and rename it to ```chr1_train.vcf.gz``` for the code to actually do anything. From there, you should just be able to run ```python imputerSE.py``` to run the main pipeline. 

As you can tell, my project is still not complete. So far, I have only been able to impute on very small test sets, and I don't have a good grasp on how I could extend the project to reasonably run on larger datasets. I will have to look into more state-of-the-art algorithms to see if I can improve performance in that regard. 

My number one priority is to write code to calculate the accuracy of my predictions, then benchmark against an industry-standard imputation tool, to see how my code performs in comparison (runtime, accuracy, etc.)

Also, one of my goals for the project was being able to predict missing SNP genotypes on my own AncestryDNA dataset. However, the SNP data from Ancestry or 23andMe is very sparse compared to the 1000Genomes training set, so there can be large gaps (multiple kb), so I would guess that the power to impute with LD is greatly reduced. I will have to decide how to approach that, since I also don't have a very accurate way to benchmark it.


