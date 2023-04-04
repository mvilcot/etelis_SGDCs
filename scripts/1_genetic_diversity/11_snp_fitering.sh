
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








#### ---- Infos ---------------------------------------------------------------------------------------------------------------------
This machine is composed of 64 on the host, and 512GB of RAM (it is possible that we access to some resources for tests purposes).

The machine is accessible via SSH to the following address:
162.38.181.160

ssh mvilcot@162.38.181.160
or
ssh mvilcot@bigmem0.mbb.univ-montp2.fr

You have a sudo access to this machine (*), so you can install all the programs you need. However, please note that you have to use it under the terms of the University (1) and renater (2) policy usage.


(*) The root password that was assigned to you is: Hfack4CX


Because of the full access given to you, it is assumed that you know some linux basics: how to use linux and how to install or remove software (although many bioinformatics software are already present (folder /mnt/bin (corresponding to /share/apps/bin on the cluster) ...)).

A large disk is mounted on /media/bigvol. Therefore, it is strongly recommended that you store your data at that location; especially in case of failure, there is high chance that only the data stored there are recoverable. However, we remind you that as for the cluster, we do not provide the backup of your data. Therefore, thank you to recover steadily your datas in the past and especially to recover BEFORE the end of the reservation. Indeed, we resettle, erase and clean the machine and the data once the booking is completed for the next booking.


For any inquiry, thank you for going through the usual interface:
http://kimura.univ-montp2.fr/aide/index.php?a=add


For any work using the service computing platform and bioinformatics Bioinformatics Biodiversity Montpellier, thank you to include this formula in your publications:

[Replace_with_your_project_name] benefited from the Montpellier Bioinformatics Biodiversity platform supported by the LabEx CeMEB, an ANR "Investissements d'avenir" program (ANR-10-LABX-04-01).

kind regards,

The HPC platform team

1] : https://www.umontpellier.fr/wp-content/uploads/2014/07/CHARTE-USAGE-SI-UMontpellier.pdf
[2] : https://www.renater.fr/telechargement,1392


