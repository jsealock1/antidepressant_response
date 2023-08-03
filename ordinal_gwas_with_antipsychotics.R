
## make file
## 1 = responder
## 2 = intermediate ad only
## 3 = intermediate + antipsychotic 
## 4 = non-responder
setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas")
library(data.table)

first_trial <- read.table("../ad_response/first_drug_trial_final_response_sd_pull_052623_2023-06-21.txt", header=T, sep="\t") 

demo = fread("/data/davis_lab/shared/phenotype_data/biovu/delivered_data/sd_wide_pull/20230607_pull/20230607_sd_pull_person_all.txt.gz")
demo$DOB = as.Date(demo$birth_datetime, "%Y-%m-%d")
dob = demo[,c("GRID","DOB")]
colnames(dob) <- c("GRID","DOB")

first_trial_dob <- merge(first_trial, dob, by="GRID")
first_trial_dob$age_at_response <- as.numeric(as.Date(first_trial_dob$last_drug_date, "%Y-%m-%d") - as.Date(first_trial_dob$DOB, "%Y-%m-%d"))/365.25
pcs <- read.table("/data/davis_lab/shared/genotype_data/biovu/processed/imputed/best_guess/MEGA/MEGA_recalled/20200330_MEGA.GRID.RACE.ETH.GEN.batch.PCs.covariates_EU.filt1.txt", header=T, sep="\t")[c(1,2,8:17)]
dat = merge(first_trial_dob, pcs, by="GRID")

dat$drug_response2 = ifelse(dat$drug_response == "responder", "1", 
						ifelse(dat$drug_response == "intermediate_antidepressants_only", "2",
							ifelse(dat$drug_response == "intermediate_multiple_antidepressant_and_antipsychotic" | dat$drug_response=="intermediate_antidepressant1_and_antipsychotic", "3", "4")))
dat2 = dat[c(1,7:19)]
dat2$GRID = paste0(dat2$GRID,"_",dat2$GRID)
write.csv(dat2, "ordinal_response_gwas_dat_4_responses_with_antipsychotics_2023-06-21.csv", col.names=T, row.names=F, quote=F)

###
# 1. get samples from vcf files
ml GCC/10.2.0 BCFtools/1.12
bcftools query -l "/data/davis_lab/shared/genotype_data/biovu/processed/imputed/dosage/MEGA/MEGA_recalled_with_AIMS/chr22.dose.all.sets.vcf.gz" > "mega_dosage_sample_ids.txt"

## 2. find overlap between ad response and vcf file samples
dat <- read.csv("ordinal_response_gwas_dat_4_responses_with_antipsychotics_2023-06-21.csv", header=T)

vcf <- read.table("mega_dosage_sample_ids.txt", header=F)
nrow(dat)
# [1] 30152
nrow(dat[(dat$GRID %in% vcf$V1),])
# [1] 30152

### all ad response samples are in vcf file

samples = dat[c(1)]
write.table(samples, "ad_response_samples.txt", col.names=F, quote=F, row.names=F)


#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=250GB
#SBATCH --array=1-22
#SBATCH --cpus-per-task 3
#SBATCH --time=96:00:00
#SBATCH --output=/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/files/get_files.txt
#SBATCH --mail-user=jsealock@broadinstitute.org
#SBATCH --mail-type=ALL
#SBATCH --job-name="mega_dosage_snp_qc"
#SBATCH --account=davis_lab

cd /data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/files/

ml GCC/10.2.0 BCFtools/1.12

bcftools view -T mega_qced_snps2.txt --samples-file ad_response_samples.txt -Oz -o ad_response_samples_dosage_chr${SLURM_ARRAY_TASK_ID}.filtered.vcf.gz /data/davis_lab/shared/genotype_data/biovu/processed/imputed/dosage/MEGA/MEGA_recalled_with_AIMS/chr${SLURM_ARRAY_TASK_ID}.dose.all.sets.vcf.gz

sbatch /data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/mega_gwas.slurm

#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=250GB
#SBATCH --array=1,2,3,4,11,12
#SBATCH --cpus-per-task 3
#SBATCH --time=200:00:00
#SBATCH --output=/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/mega_ordinal_gwas_LRT_output.txt
#SBATCH --mail-user=jsealock@broadinstitute.org
#SBATCH --mail-type=ALL
#SBATCH --job-name="mega_ord_lrt_gwas"
#SBATCH --account=davis_lab

module load Julia/1.7.2
cd /data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/

julia mega_gwas_lrt_script.jl ${SLURM_ARRAY_TASK_ID}




### GWAS script
chr_num = ARGS[1]

print(chr_num)

vcffile = "/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/files/ad_response_samples_dosage_chr"*chr_num*".filtered"

print(vcffile)

covfile = "/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/ordinal_response_gwas_dat_4_responses_with_antipsychotics_2023-06-21.csv"

pvaloutputfile = "/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/"*chr_num*"ordinal_response_gwas_dat_4_responses_with_antipsychotics"


using BenchmarkTools, CSV, Glob, SnpArrays, OrdinalGWAS, DataFrames

nm = ordinalgwas(@formula(drug_response2 ~ age_at_response + GENDER + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10), covfile, nothing)

ordinalgwas(nm, vcffile, test=:LRT, geneticformat = "VCF", vcftype = :DS, pvalfile=pvaloutputfile)


###################################################################################################
# cat *ordinal_response_gwas_dat_4_responses_with_antipsychotics > ordinal_response_gwas_dat_4_responses_with_antipsychotics_results.txt

library(data.table)
library(tidyr)
library(qqman)
# setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/")

chr1 = read.csv("1ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr2 = read.csv("2ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr3 = read.csv("3ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr4 = read.csv("4ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr5 = read.csv("5ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr6 = read.csv("6ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr7 = read.csv("7ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr8 = read.csv("8ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr9 = read.csv("9ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr10 = read.csv("10ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr11 = read.csv("11ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr12 = read.csv("12ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr13 = read.csv("13ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr14 = read.csv("14ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr15 = read.csv("15ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr16 = read.csv("16ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr17 = read.csv("17ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr18 = read.csv("18ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr19 = read.csv("19ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr20 = read.csv("20ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr21 = read.csv("21ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)
chr22 = read.csv("22ordinal_response_gwas_dat_4_responses_with_antipsychotics", header=T)

gwas = rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22)

# gwas <- fread("ordinal_response_gwas_dat_4_responses_with_antipsychotics_all.txt", header=T)
gwas$pval <- as.numeric(gwas$pval)
gwas$pos <- as.numeric(gwas$pos)
gwas$chr <- as.numeric(gwas$chr)
# gwas_complete <- gwas[complete.cases(gwas),]

gwas2 <- separate(gwas, snpid, c("chr2","pos2","ref","alt"), sep=":")
gwas2$loc = paste0(gwas2$chr2, ":", gwas2$pos2)
gwas3 = gwas2[c(3:10)]
colnames(gwas3)[1] = "chr"
colnames(gwas3)[2] = "pos"


## add rsids
bim = read.table("/data/davis_lab/sealockj/projects/scripts/mega_geno_data/20200330_biallelic_mega_recalled.chr1-22.grid.EU.filt1_IBD_filt_0.2.bim", header=F)
bim$loc = paste0(bim$V1, ":", bim$V4)

# gwas_loc_dedup_hardcall_snps = gwas_loc_dedup[(gwas_loc_dedup$snp %in% bim$V2),]
gwas_snps = merge(gwas3, bim, by="loc")
gwas_snps2 = gwas_snps[c(10,2:8)]

colnames(gwas_snps2)[1] = "rsid"
write.table(gwas_snps2, "ad_response_mega_ordinal_gwas_hardcalled_qced_snps.txt",  col.names=T, row.names=F, quote=F, sep="\t")


## 

dat = read.table("ad_response_mega_ordinal_gwas_hardcalled_qced_snps.txt", header=T, sep="\t")

dat = gwas_snps2

png("ad_response_lrt_gwas_mega_dosage_qqplot.png")
p = qq(dat$pval)
print(p)
dev.off()

png("ad_response_lrt_gwas_mega_dosage_manhattan.png")
p = manhattan(dat, chr="chr", bp="pos", p="pval", snp="rsid")
print(p)
dev.off()


