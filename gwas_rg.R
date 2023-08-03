# load results
library(data.table)
library(tidyr)
library(qqman)
setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/")

gwas <- fread("ordinal_response_gwas_dat_4_responses_with_antipsychotics_results.txt", header=T)
gwas$pval <- as.numeric(gwas$pval)
gwas$pos <- as.numeric(gwas$pos)
gwas$chr <- as.numeric(gwas$chr)
gwas_complete <- gwas[complete.cases(gwas),]

gwas2 <- separate(gwas_complete, snpid, c("chr2","pos2","ref","alt"), sep=":")
gwas2$chr2 <- NULL
gwas2$pos2 <- NULL

bim = read.table("/data/davis_lab/sealockj/projects/scripts/mega_geno_data/20200330_biallelic_mega_recalled.chr1-22.grid.EU.filt1_IBD_filt_0.2.bim", header=F)

gwas2$loc = paste0(gwas2$chr, ":", gwas2$pos)
bim$loc = paste0(bim$V1, ":", bim$V4)
bim2 = bim[c(2,7)]
gwas3 = merge(gwas2, bim2, by="loc")


nrow(bim)
# [1] 36,604,696

nrow(gwas)
# [1] 12,308,366

nrow(gwas3)
# [1] 9,972,110

loc = read.table("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/pgs_assoc/gwas/openmendel/dosage/mega/hg19_avsnp147.txt.gz", header=F, sep="\t")
loc$loc <- paste0(loc$V1,":",loc$V2)
loc2 <- loc[c(6,7)]
gwas4 = merge(gwas2, loc2, by="loc")
# [1] 12,940,260

colnames(gwas4)[9] <- 'snp'
gwas4$N = 30152
write.table(gwas4, '/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/ad_response_mega_ordinal_gwas_all_snps.txt', row.names=F, col.names=T, quote=F, sep="\t")


######
### calculate rg
#### rg
ad_gwas_file="/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/ad_response_mega_ordinal_gwas_all_snps.txt"
ld_path="/data/davis_lab/sealockj/projects/cad_labwas_ms/sdwide_current/hdl_method_rg/HDL/UKB_imputed_SVD_eigen99_extraction"
working_dir="/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/gwas/genetic_correlation/"


## wrangle sumstats 
# ad response gwas
Rscript /data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/pgs_assoc/gwas/saige/hdl_genetic_cor/HDL/HDL.data.wrangling.R \
gwas.file=$ad_gwas_file \
LD.path=$ld_path \
SNP=snp A1=alt A2=ref N=N b=effect se=stder \
output.file="$working_dir/ad_response" \
log.file="$working_dir/ad_response"


###################################################################################################

## calculate rg 
psych_dir="/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/pgs_assoc/gwas/saige/hdl_genetic_cor/"

## depression 
Rscript /data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sdwide/gwas/results/HDL/HDL.run.R \
gwas1.df="$working_dir/ad_response.hdl.rds" \
gwas2.df="$psych_dir/depression.hdl.rds" \
LD.path=$ld_path \
output.file="$working_dir/depression_vs_ad_response"

## anxiety 
Rscript /data/davis_lab/sealockj/projects/cad_labwas_ms/sdwide_current/hdl_method_rg/HDL/HDL.run.R \
gwas1.df="$working_dir/ad_response.hdl.rds" \
gwas2.df="$psych_dir/anx.hdl.rds" \
LD.path=$ld_path \
output.file="$working_dir/anx_vs_ad_response"

#bip
Rscript /data/davis_lab/sealockj/projects/cad_labwas_ms/sdwide_current/hdl_method_rg/HDL/HDL.run.R \
gwas1.df="$working_dir/ad_response.hdl.rds" \
gwas2.df="$psych_dir/bip.hdl.rds" \
LD.path=$ld_path \
output.file="$working_dir/bip_vs_ad_response"

#scz
Rscript /data/davis_lab/sealockj/projects/cad_labwas_ms/sdwide_current/hdl_method_rg/HDL/HDL.run.R \
gwas1.df="$working_dir/ad_response.hdl.rds" \
gwas2.df="$psych_dir/scz.hdl.rds" \
LD.path=$ld_path \
output.file="$working_dir/scz_vs_ad_response"

## cross disorders
Rscript /data/davis_lab/sealockj/projects/cad_labwas_ms/sdwide_current/hdl_method_rg/HDL/HDL.run.R \
gwas1.df="$working_dir/ad_response.hdl.rds" \
gwas2.df="$psych_dir/cross_disorders.hdl.rds" \
LD.path=$ld_path \
output.file="$working_dir/cross_disorders_vs_ad_response"

###################################################################################################

library(ggplot2)
library(devtools)
library(ggsci)

dat <- read.csv("rg_results.csv", header=T)

dat$Lower.CI <- as.numeric(dat$rg - (1.96*dat$SE))
dat$Upper.CI <- as.numeric(dat$rg + (1.96*dat$SE))

dat3 <- subset(dat, Trait!="Anxiety")
dat3$Trait <- factor(dat3$Trait, levels=c("Cross Disorders","Schizophrenia","Bipolar","Depression"))

pdf("ad_response_rg_mega_ordinal.pdf", width=8, height=5)
p <- ggplot(data=dat3, aes(x=Trait, y=rg, ymin=Lower.CI, ymax=Upper.CI, fill=Trait, colour=Trait)) + 
    geom_pointrange(position=position_dodge(width=0.01), size=1.5) + 
    geom_hline(yintercept=0, lty=2) + 
    coord_flip() + 
    ggtitle(" ") +
    theme_classic() +
    xlab(" ") +
    theme(axis.text.x=element_text(size=16), axis.title.x=element_text(size=16)) +
    theme(legend.title = element_text(size = 16), legend.text = element_text(size=16), axis.text.y=element_text(size=16)) +
    ylab("Genetic Correlation") +
     scale_color_manual(values=c("Depression" = "#F39B7FFF", "Bipolar" = "#00A087FF", "Schizophrenia" = "#4DBBD5FF", "Cross Disorders" = "#DC0000FF")) +
    guides(colour=F, fill=F)
print(p)
dev.off()




