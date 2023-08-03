### association between phq and ad response
library(data.table)
library(lme4)
library(ggplot2)
library(dplyr)
library(ggsci)
date = Sys.Date()

setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/phq_data/")

## load data
phq8 = read.table("phq8_clean_data_sd_wide_pull_20230607.txt", header=T, sep="\t")
phq2 = read.table("phq2_clean_data_sd_wide_pull_20230607.txt", header=T, sep="\t")

ad_response <- read.table("../ad_response/antidepressant_response_output_with_antipsychotics_052623_sd_pull_062123.txt", header=T, sep="\t")

## try new demo data
demo = fread("/data/davis_lab/shared/phenotype_data/biovu/delivered_data/sd_wide_pull/20230607_pull/20230607_sd_pull_person_all.txt.gz")
demo$DOB = as.Date(demo$birth_datetime, "%Y-%m-%d")
demo = demo[,c("GRID","gender_source_value","DOB","race_source_value")]
colnames(demo) <- c("GRID","GENDER","DOB","RACE")

## organize data frames
# ad_response2 = ad_response[c(1:4,8,19)]
ad_response2 = ad_response[,c("GRID","Short_Name","MED_DATE","drug_response")]

## find n individuals with response and any phq 
ad_phq8 <- merge(ad_response2, phq8, by="GRID")
nrow(ad_phq8[!duplicated(ad_phq8$GRID),])
#[1] 15319

ad_phq2 <- merge(ad_response2, phq2, by="GRID")
nrow(ad_phq2[!duplicated(ad_phq2$GRID),])
# [1] 25964

min(as.Date(phq2$date, "%Y-%m-%d"))
# [1] "2018-08-16"
min(as.Date(phq8$date, "%Y-%m-%d"))
# [1] "2017-08-10"


## find phq closest to ad response date
create_long_model_data_all <- function(input=input, phq=phq, time=time){
    phq$phq_date <- as.Date(phq$date, "%Y-%m-%d")
    phq_numbered <- phq %>%
        group_by(GRID) %>%
        arrange(phq_date) %>%
        mutate(observation=1:n())
    ad_response_and_phq <- merge(input, phq_numbered, by="GRID")
    ad_response_and_phq$time_diff <- as.numeric(as.Date(ad_response_and_phq$phq_date, "%Y-%m-%d") - as.Date(ad_response_and_phq$MED_DATE, "%Y-%m-%d"))
    ## for each individual and PHQ observations, find closest med date
    closest_date <- setDT(ad_response_and_phq)[, .SD[which.min(abs(time_diff))], by=c("GRID","observation")]
    closest_date_time <- subset(closest_date, abs(time_diff)<=as.numeric(time))
    out <- merge(closest_date_time, demo, by="GRID")
    return(out)
}

## run long models
## run longitudinal model PHQ ~ AD response within class 
### use PHQ that come within 2 weeks of AD respone variable 

long_model <- function(dat=dat, med=med, version=version){
    dat$drug_response = as.character(dat$drug_response)
    dat$age_at_event <- as.numeric(as.Date(dat$phq_date, "%Y-%m-%d") - as.Date(dat$DOB, "%Y-%m-%d"))/365.25
    dat$age_at_event <- scale(dat$age_at_event)
    dat$answer <- scale(dat$answer)
    if(version==4){
        dat$drug_response <- ifelse(dat$drug_response == "intermediate_multiple_antidepressant_and_antipsychotic" | dat$drug_response ==  "intermediate_antidepressant1_and_antipsychotic", "intermediate_with_antipsychotic", dat$drug_response)
        dat$drug_response <- ifelse(dat$drug_response == "nonresponder_no_antipsychotics" | dat$drug_response == "nonresponder_with_antipsychotics", "nonresponder", dat$drug_response)
        dat$drug_response <- factor(dat$drug_response, levels=c("responder","intermediate_antidepressants_only","intermediate_with_antipsychotic", "nonresponder"))
    }
    if(version==6){
    dat$drug_response <- factor(dat$drug_response, levels=c("responder","intermediate_antidepressants_only","intermediate_antidepressant1_and_antipsychotic", "intermediate_multiple_antidepressant_and_antipsychotic", "nonresponder_no_antipsychotics", "nonresponder_with_antipsychotics"))
    }
    tvc.model <- lmer(answer ~ drug_response + RACE + GENDER + age_at_event + 
        (1|GRID), data=dat, REML=F)
    null.model <- lmer(answer ~ RACE + GENDER + age_at_event + (1|GRID), data=dat, REML=F)
    anova <- anova(tvc.model, null.model)
    pvalue = anova[2,8]
    tvc.effect = coef(summary(tvc.model))[2,1]
    tvc.se = coef(summary(tvc.model))[2,2]
    n_obs = nrow(dat)
    n_inds = nrow(dat[!duplicated(dat$GRID),])
    out <- data.frame(med, version, pvalue, tvc.effect, tvc.se, n_obs, n_inds)
}

# ## create models 

# all classes
phq8_all_classes <- create_long_model_data_all(input=ad_response2, phq=phq8, time=14)
# find individuals with multiple instances in the dataset
phq8_all_classes = phq8_all_classes %>% group_by(GRID) %>% mutate(observation=1:n())
phq8_all_classes_mult = subset(phq8_all_classes, observation>1)

phq8_all_classes_mult2 = phq8_all_classes[(phq8_all_classes$GRID %in% phq8_all_classes_mult$GRID),]

all_classes_long_phq8_6_responses <- long_model(dat=phq8_all_classes_mult2, med="PHQ8", version=6)
all_classes_long_phq8_4_responses <- long_model(dat=phq8_all_classes_mult2, med="PHQ8", version=4)
all_classes_long_phq8 = rbind(all_classes_long_phq8_4_responses, all_classes_long_phq8_6_responses)




## phq2 

phq2_all_classes <- create_long_model_data_all(input=ad_response2, phq=phq2, time=14)
# find individuals with multiple instances in the dataset
phq2_all_classes = phq2_all_classes %>% group_by(GRID) %>% mutate(observation=1:n())
phq2_all_classes_mult = subset(phq2_all_classes, observation>1)

phq2_all_classes_mult2 = phq2_all_classes[(phq2_all_classes$GRID %in% phq2_all_classes_mult$GRID),]

all_classes_long_phq2_6_responses <- long_model(dat=phq2_all_classes_mult2, med="PHQ2", version=6)
all_classes_long_phq2_4_responses <- long_model(dat=phq2_all_classes_mult2, med="PHQ2", version=4)
all_classes_long_phq2 = rbind(all_classes_long_phq2_4_responses, all_classes_long_phq2_6_responses)

# save
out = rbind(all_classes_long_phq8, all_classes_long_phq2)
date = Sys.Date()
write.csv(out, paste0("phq_ad_response_with_antipsychotics_multiple_observations_longitudinal_results_all_samples_",date,".csv"), row.names=F, quote=F)


# r2
write.table(phq8_all_classes_mult2, "ad_response_phq8_model_dat.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(phq2_all_classes_mult2, "ad_response_phq2_model_dat.txt", col.names=T, row.names=F, quote=F, sep="\t")

dat$drug_response = as.character(dat$drug_response)
dat$age_at_event <- as.numeric(as.Date(dat$phq_date, "%Y-%m-%d") - as.Date(dat$DOB, "%Y-%m-%d"))/365.25
dat$age_at_event <- scale(dat$age_at_event)
dat$answer <- scale(dat$answer)
    dat$drug_response <- ifelse(dat$drug_response == "intermediate_multiple_antidepressant_and_antipsychotic" | dat$drug_response ==  "intermediate_antidepressant1_and_antipsychotic", "intermediate_with_antipsychotic", dat$drug_response)
    dat$drug_response <- ifelse(dat$drug_response == "nonresponder_no_antipsychotics" | dat$drug_response == "nonresponder_with_antipsychotics", "nonresponder", dat$drug_response)
    dat$drug_response <- factor(dat$drug_response, levels=c("responder","intermediate_antidepressants_only","intermediate_with_antipsychotic", "nonresponder"))

tvc.model <- lmer(answer ~ drug_response + RACE + GENDER + age_at_event + (1|GRID), data=dat, REML=F)
null.model <- lmer(answer ~ RACE + GENDER + age_at_event + (1|GRID), data=dat, REML=F)

library(MuMIn)
r.squaredGLMM(tvc.model)
r.squaredGLMM(null.model)


#####################
# subset to depression cases 
phecode = readRDS("/data/davis_lab/shared/phenotype_data/biovu/phecode_tables/medical_home_phecode_table_20210806_pull_remake_111822_2_distinct_dates_no_exclusions.Rds")
# 296.22 mdd
# 296.2 dep 
# 304 - adjustment rxn 
# 300.4 - dysthymic 

dep = phecode[c(1,487,488,512,499)]
colnames(dep) <- c("GRID","dep","mdd","adj_rxn","dysthymic")
dep_cases = subset(dep, dep=="TRUE" | mdd=="TRUE" | adj_rxn=="TRUE" | dysthymic=="TRUE")

# phq8 
phq8_all_classes_mult_dep = phq8_all_classes_mult2[(phq8_all_classes_mult2$GRID %in% dep_cases$GRID),]
all_classes_long_phq8_6_responses_dep <- long_model(dat=phq8_all_classes_mult_dep, med="PHQ8", version=6)
all_classes_long_phq8_4_responses_dep <- long_model(dat=phq8_all_classes_mult_dep, med="PHQ8", version=4)
all_classes_long_phq8_dep = rbind(all_classes_long_phq8_4_responses_dep, all_classes_long_phq8_6_responses_dep)

# phq2
phq2_all_classes_mult_dep = phq2_all_classes_mult2[(phq2_all_classes_mult2$GRID %in% dep_cases$GRID),]
all_classes_long_phq2_6_responses_dep <- long_model(dat=phq2_all_classes_mult_dep, med="PHQ2", version=6)
all_classes_long_phq2_4_responses_dep <- long_model(dat=phq2_all_classes_mult_dep, med="PHQ2", version=4)
all_classes_long_phq2_dep = rbind(all_classes_long_phq2_4_responses_dep, all_classes_long_phq2_6_responses_dep)

out = rbind(all_classes_long_phq8_dep, all_classes_long_phq2_dep)
write.csv(out, paste0("phq_ad_response_with_antipsychotics_multiple_observations_longitudinal_results_depression_cases_",date,".csv"), row.names=F, quote=F)


#####################
# plots


plot_phq = function(phq8=phq8, phq2=phq2, response_version=response_version, sample=sample){
    phq8$PHQ = "PHQ 8"
    phq2$PHQ = "PHQ 2"
    dat = rbind(phq2, phq8)
    dat$response_version <- paste0(response_version," Responses")
    if(response_version==6){
        dat$response = ifelse(dat$drug_response == "responder","Responder", 
                            ifelse(dat$drug_response == "intermediate_antidepressants_only", "Intermediate AD Only", 
                                ifelse(dat$drug_response == "intermediate_antidepressant1_and_antipsychotic", "Intermediate 1 AD + AP",
                                    ifelse(dat$drug_response == "intermediate_multiple_antidepressant_and_antipsychotic", "Intermediate >1 AD + AP", 
                                        ifelse(dat$drug_response=="nonresponder_no_antipsychotics", "Non-responder No AP", "Non-responder with AP")))))
        dat$response = factor(dat$response, levels=c("Responder", "Intermediate AD Only", "Intermediate 1 AD + AP", "Intermediate >1 AD + AP", "Non-responder No AP", "Non-responder with AP"))
    }
    if(response_version==4){
        dat$response = ifelse(dat$drug_response == "responder","Responder", 
                            ifelse(dat$drug_response == "intermediate_antidepressants_only", "Intermediate AD Only", 
                                ifelse(dat$drug_response == "intermediate_antidepressant1_and_antipsychotic", "Intermediate AD + AP",
                                    ifelse(dat$drug_response == "intermediate_multiple_antidepressant_and_antipsychotic", "Intermediate AD + AP", "Non-responder"))))
        dat$response = factor(dat$response, levels=c("Responder", "Intermediate AD Only", "Intermediate AD + AP", "Non-responder"))
    }

    pdf(paste0("phq_ad_response_longitudinal_results_",response_version,"_responses","_",sample,"_",date,".pdf"), width=10, height=8)
    p8 <- ggplot(dat, aes(x=response, y=answer, fill=response)) +
        geom_boxplot() +
        theme_bw() +
        facet_wrap(~ PHQ) +
        theme(axis.text.x = element_text(angle = 45, hjust=1), plot.margin = margin(10, 10, 10, 100)) +
        ylab("PHQ Score") +
        scale_fill_npg() +
        theme(axis.text.x=element_blank()) +
        xlab("Antidepressant Response") +
        guides(fill=guide_legend(title="Antidepressant Response"))
    print(p8)
    dev.off()
}
plot_phq(phq8=phq8_all_classes_mult, phq2=phq2_all_classes_mult2, response_version=6, sample="all")
plot_phq(phq8=phq8_all_classes_mult, phq2=phq2_all_classes_mult2, response_version=4, sample="all")

plot_phq(phq8=phq8_all_classes_mult_dep, phq2=phq2_all_classes_mult_dep, response_version=6, sample="dep_cases")
plot_phq(phq8=phq8_all_classes_mult_dep, phq2=phq2_all_classes_mult_dep, response_version=4, sample="dep_cases")


plot_manuscript = function(dat=dat, phq_version=phq_version){
    dat$PHQ = paste0("PHQ ", phq_version)
    dat$response = ifelse(dat$drug_response == "responder","Responder", 
                            ifelse(dat$drug_response == "intermediate_antidepressants_only", "Intermediate AD Only", 
                                ifelse(dat$drug_response == "intermediate_antidepressant1_and_antipsychotic", "Intermediate AD + AP",
                                    ifelse(dat$drug_response == "intermediate_multiple_antidepressant_and_antipsychotic", "Intermediate AD + AP", "Non-responder"))))

    dat$response = factor(dat$response, levels=c("Responder", "Intermediate AD Only", "Intermediate AD + AP", "Non-responder"))

    pdf(paste0("phq",phq_version,"_ad_response_longitudinal_results_4_responses","_",date,".pdf"), width=8, height=6)
    p8 <- ggplot(dat, aes(x=response, y=answer, fill=response)) +
        geom_boxplot() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust=1), plot.margin = margin(10, 10, 10, 100)) +
        ylab("PHQ Score") +
        scale_fill_npg() +
        theme(axis.text.x=element_blank()) +
        xlab("Antidepressant Response") +
        guides(fill=guide_legend(title="Antidepressant Response")) +
        facet_wrap(~sample)
    print(p8)
    dev.off()
}


phq8_all_classes_mult2$sample = "All"
phq8_all_classes_mult_dep$sample = "Depression Cases"
phq8_dat = rbind(phq8_all_classes_mult2, phq8_all_classes_mult_dep)
plot_manuscript(dat=phq8_dat, phq_version="8")

phq2_all_classes_mult2$sample = "All"
phq2_all_classes_mult_dep$sample = "Depression Cases"
phq2_dat = rbind(phq2_all_classes_mult2, phq2_all_classes_mult_dep)
plot_manuscript(dat=phq2_dat, phq_version="2")


plot_manuscript2 = function(dat=dat, phq_version=phq_version){
    dat$PHQ = paste0("PHQ ", phq_version)
    dat$response = ifelse(dat$drug_response == "responder","Responder", 
                        ifelse(dat$drug_response == "intermediate_antidepressants_only", "Intermediate AD Only", 
                            ifelse(dat$drug_response == "intermediate_antidepressant1_and_antipsychotic", "Intermediate 1 AD + AP",
                                ifelse(dat$drug_response == "intermediate_multiple_antidepressant_and_antipsychotic", "Intermediate >1 AD + AP", 
                                    ifelse(dat$drug_response=="nonresponder_no_antipsychotics", "Non-responder No AP", "Non-responder with AP")))))
    dat$response = factor(dat$response, levels=c("Responder", "Intermediate AD Only", "Intermediate 1 AD + AP", "Intermediate >1 AD + AP", "Non-responder No AP", "Non-responder with AP"))

    pdf(paste0("phq",phq_version,"_ad_response_longitudinal_results_6_responses","_",date,".pdf"), width=8, height=6)
    p8 <- ggplot(dat, aes(x=response, y=answer, fill=response)) +
        geom_boxplot() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust=1), plot.margin = margin(10, 10, 10, 100)) +
        ylab("PHQ Score") +
        scale_fill_npg() +
        theme(axis.text.x=element_blank()) +
        xlab("Antidepressant Response") +
        guides(fill=guide_legend(title="Antidepressant Response")) +
        facet_wrap(~sample)
    print(p8)
    dev.off()
}

plot_manuscript2(dat=phq8_dat, phq_version="8")
plot_manuscript2(dat=phq2_dat, phq_version="2")

