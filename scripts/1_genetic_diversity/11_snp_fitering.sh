
#### ---- Script ---------------------------------------------------------------------------------------------------------------------


ssh mvilcot@162.38.181.160
ssh mvilcot@bigmem0.mbb.univ-montp2.fr
Hfack4CX

# ~/R/x86_64-pc-linux-gnu-library/3.6
# module load modulefiles/R/4.0.5


### 0. Raw file
# 105144 sites, 369 samples
Data/DartSeq/Report_DEtel22-6705_SNP.vcf


## 1. Individuals with too much missings
# 105144 SNPs, 364 samples
# Asses individual levels of missing data
vcftools --vcf Data/DartSeq/Report_DEtel22-6705_SNP.vcf --missing-indv
# Print percentage of missing by ind
cat out.imiss 
# Create a list of inidivduals with more than 50% missings
awk '($5 > 0.5)' out.imiss | cut -f1 > Intermediate/lowDP.indv
# Remove those individuals from vcf 
vcftools --vcf Data/DartSeq/Report_DEtel22-6705_SNP.vcf --remove Intermediate/lowDP.indv --recode --recode-INFO-all --out Intermediate/Report_DEtel22-6705_SNP_missind


## 2. Missings by site
# 87641 SNPs 
vcftools --vcf Intermediate/Report_DEtel22-6705_SNP_missind.recode.vcf --max-missing 0.70 --recode --recode-INFO-all --out Intermediate/Report_DEtel22-6705_SNP_missind_callrate0.70


## 3. Minimum allele frequency 
# 21948 SNPs
vcftools --vcf Intermediate/Report_DEtel22-6705_SNP_missind_callrate0.70.recode.vcf --maf 0.05 --recode --recode-INFO-all --out Intermediate/Report_DEtel22-6705_SNP_missind_callrate0.70_maf0.05


## Heterozygosity
vcftools --vcf Intermediate/Report_DEtel22-6705_SNP_missind_callrate0.70.recode.vcf --geno-r2



