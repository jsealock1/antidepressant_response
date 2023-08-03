
setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/pgs_assoc/")

library(data.table)
library(MASS)

## associate with response 
first_trial = read.table("../ad_response/first_drug_trial_final_response_sd_pull_052623_2023-06-21.txt", header=T)

demo = fread("/data/davis_lab/shared/phenotype_data/biovu/delivered_data/sd_wide_pull/20230607_pull/20230607_sd_pull_person_all.txt.gz")
demo$DOB = as.Date(demo$birth_datetime, "%Y-%m-%d")
dob = demo[,c("GRID","DOB")]
colnames(dob) <- c("GRID","DOB")

first_trial_dob <- merge(first_trial, dob, by.x="GRID")
first_trial_dob$age_at_response <- as.numeric(as.Date(first_trial_dob$last_drug_date, "%Y-%m-%d") - as.Date(first_trial_dob$DOB, "%Y-%m-%d"))/365.25

date = Sys.Date()
ordinal_regression = function(dat=dat, name=name, n_responses=n_responses){
    dat$drug_response = as.character(dat$drug_response)
    # dat$nonremission <- scale(dat$nonremission)
    if(n_responses == 6){
        dat$drug_response <- factor(dat$drug_response, levels=c("responder", "intermediate_antidepressants_only", "intermediate_antidepressant1_and_antipsychotic", "intermediate_multiple_antidepressant_and_antipsychotic","nonresponder_no_antipsychotics","nonresponder_with_antipsychotics"))
    }
    if(n_responses == 4){
        dat$drug_response <- ifelse(dat$drug_response == "intermediate_multiple_antidepressant_and_antipsychotic" | dat$drug_response ==  "intermediate_antidepressant1_and_antipsychotic", "intermediate_with_antipsychotic", dat$drug_response)
        dat$drug_response <- ifelse(dat$drug_response == "nonresponder_no_antipsychotics" | dat$drug_response == "nonresponder_with_antipsychotics", "nonresponder", dat$drug_response)

        dat$drug_response <- factor(dat$drug_response, levels=c("responder", "intermediate_antidepressants_only","intermediate_with_antipsychotic", "nonresponder"))
    }

    dat$anx <- scale(dat$anx)
    dat$bip <- scale(dat$bip)
    dat$cross_disorders <- scale(dat$cross_disorders)
    dat$scz <- scale(dat$scz)
    dat$dep <- scale(dat$dep)

    print("run full model")
    output <- polr(drug_response ~ anx + bip + dep + scz + age_at_response + GENDER + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat, Hess=TRUE)
    ctable <- coef(summary(output))
    p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
    ctable2 <- cbind(ctable, "p value" = p)
    write.table(ctable2, paste0(name, "_stats_summary_first_response_vs_psych_pgs_",n_responses, "_responses_",date,".txt"), row.names=T, col.names=T, sep="\t", quote=F)

    print("run cd model")
    cd_output <- polr(drug_response ~ cross_disorders + age_at_response + GENDER + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat, Hess=TRUE)
    cd_ctable <- coef(summary(cd_output))
    cd_p <- pnorm(abs(cd_ctable[, "t value"]), lower.tail = FALSE) * 2
    cd_ctable2 <- cbind(cd_ctable, "p value" = cd_p)
    print("writing cd results")
    write.table(cd_ctable2, paste0(name, "_stats_summary_first_response_vs_cross_disorders_psych_pgs_",n_responses, "_responses_",date,".txt"), row.names=T, col.names=T, sep="\t", quote=F)
}

# EUR
# nonresmission_pgs <- read.table("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/pgs_assoc/antidepressant_response_pgs/nonresmission/prs_cs/ad_nonremission_prs_cs_results.profile", header=T)[c(2,6)]
# colnames(nonresmission_pgs)[2] <- "nonremission"
psych_pgs <- read.table("/data/davis_lab/sealockj/psych_pgs/psychiatric_polygenic_scores_european_ancestry_05-06-2022.txt", header=T, sep="\t")

# pgs <- merge(psych_pgs, nonresmission_pgs, by="IID")
pcs <- read.table("/data/davis_lab/shared/genotype_data/biovu/processed/imputed/best_guess/MEGA/MEGA_recalled/20200330_MEGA.GRID.RACE.ETH.GEN.batch.PCs.covariates_EU.filt1.txt", header=T, sep="\t")[c(1,2,8:17)]
pgs_pcs <- merge(psych_pgs, pcs, by.x="IID", by.y="GRID")

eur_dat <- merge(first_trial_dob, pgs_pcs, by.y="IID", by.x="GRID")
nrow(eur_dat)
# [1] 27051

ordinal_regression(dat=eur_dat, name="eur_ordinal_regression", n_responses=6)
ordinal_regression(dat=eur_dat, name="eur_ordinal_regression", n_responses=4)



#### AFR 

pcs <- read.table("/data/davis_lab/shared/genotype_data/biovu/processed/imputed/best_guess/MEGA/MEGA_recalled/20200515_MEGA.GRID.RACE.ETH.GEN.batch.PCs.covariates_AA.filt1.txt", header=T, sep="\t")[c(1,2,8:17)]
nonresmission_pgs <- read.table("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/pgs_assoc/antidepressant_response_pgs/nonresmission/afr/prs_cs/ad_nonremission_afr_prs_cs_results.profile", header=T)[c(2,6)]
psych_pgs <- read.table("/data/davis_lab/sealockj/psych_pgs/psychiatric_polygenic_scores_african_ancestry_05-06-2022.txt", header=T, sep="\t")
colnames(nonresmission_pgs)[2] <- "nonremission"
pgs <- merge(psych_pgs, nonresmission_pgs, by="IID")
pgs_pcs <- merge(pgs, pcs, by.x="IID", by.y="GRID")

afr_dat <- merge(first_trial_dob, pgs_pcs, by.y="IID", by.x="GRID")
nrow(afr_dat)
# [1] 3214

afr_dat$drug_response = as.character(afr_dat$drug_response)
ordinal_regression(dat=afr_dat, name="afr_ordinal_regression", n_responses=6)
ordinal_regression(dat=afr_dat, name="afr_ordinal_regression", n_responses=4)

################################
## control for depression dx  ##
################################

phecode = readRDS("/data/davis_lab/shared/phenotype_data/biovu/phecode_tables/medical_home_phecode_table_20210806_pull_remake_111822_2_distinct_dates_no_exclusions.Rds")
# 296.22 mdd
# 296.2 dep 
# 304 - adjustment rxn 
# 300.4 - dysthymic 

dep = phecode[c(1,487,488,512,499)]
colnames(dep) <- c("GRID","dep","mdd","adj_rxn","dysthymic")
dep_cases = subset(dep, dep=="TRUE" | mdd=="TRUE" | adj_rxn=="TRUE" | dysthymic=="TRUE")
dep_controls = subset(dep, dep=="FALSE" & mdd=="FALSE" & adj_rxn=="FALSE" & dysthymic=="FALSE")
dep_cases$depression = "TRUE"
dep_controls$depression = "FALSE"
dep = rbind(dep_cases, dep_controls)
dep = dep[c(1,6)]

eur_dat2 = merge(eur_dat, dep, by="GRID")
afr_dat2 = merge(afr_dat, dep, by="GRID")
nrow(eur_dat2)
# [1] 19814
nrow(afr_dat2)
# [1] 2444

ordinal_regression_dep_cov = function(dat=dat, name=name, n_responses=n_responses){
    dat$drug_response = as.character(dat$drug_response)
    # dat$nonremission <- scale(dat$nonremission)
    if(n_responses == 6){
        dat$drug_response <- factor(dat$drug_response, levels=c("responder", "intermediate_antidepressants_only", "intermediate_antidepressant1_and_antipsychotic", "intermediate_multiple_antidepressant_and_antipsychotic","nonresponder_no_antipsychotics","nonresponder_with_antipsychotics"))
    }
    if(n_responses == 4){
        dat$drug_response <- ifelse(dat$drug_response == "intermediate_multiple_antidepressant_and_antipsychotic" | dat$drug_response ==  "intermediate_antidepressant1_and_antipsychotic", "intermediate_with_antipsychotic", dat$drug_response)
        dat$drug_response <- ifelse(dat$drug_response == "nonresponder_no_antipsychotics" | dat$drug_response == "nonresponder_with_antipsychotics", "nonresponder", dat$drug_response)

        dat$drug_response <- factor(dat$drug_response, levels=c("responder", "intermediate_antidepressants_only","intermediate_with_antipsychotic", "nonresponder"))
    }

    dat$anx <- scale(dat$anx)
    dat$bip <- scale(dat$bip)
    dat$cross_disorders <- scale(dat$cross_disorders)
    dat$scz <- scale(dat$scz)
    dat$dep <- scale(dat$dep)

    print("run full model")
    output <- polr(drug_response ~ anx + bip + dep + scz + age_at_response + depression + GENDER + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat, Hess=TRUE)
    ctable <- coef(summary(output))
    p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
    ctable2 <- cbind(ctable, "p value" = p)
    write.table(ctable2, paste0(name, "_stats_summary_first_response_vs_psych_pgs_",n_responses, "_responses_depression_cov_",date,".txt"), row.names=T, col.names=T, sep="\t", quote=F)

    print("run cd model")
    cd_output <- polr(drug_response ~ cross_disorders + age_at_response + GENDER + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat, Hess=TRUE)
    cd_ctable <- coef(summary(cd_output))
    cd_p <- pnorm(abs(cd_ctable[, "t value"]), lower.tail = FALSE) * 2
    cd_ctable2 <- cbind(cd_ctable, "p value" = cd_p)
    print("writing cd results")
    write.table(cd_ctable2, paste0(name, "_stats_summary_first_response_vs_cross_disorders_psych_pgs_",n_responses, "_responses_depression_cov",date,".txt"), row.names=T, col.names=T, sep="\t", quote=F)
}


ordinal_regression_dep_cov(dat=eur_dat2, name="eur_ordinal_regression_dep_cov", n_responses=4)
ordinal_regression_dep_cov(dat=afr_dat2, name="afr_ordinal_regression_dep_cov", n_responses=4)
