setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607")

library(data.table)

response = read.table("ad_response/antidepressant_response_output_with_antipsychotics_052623_sd_pull_062123.txt", header=T)
first_trial <- read.table("ad_response/first_drug_trial_final_response_sd_pull_052623_2023-06-21.txt", header=T, sep="\t") 
first_trial$Short_Name = tolower(first_trial$Short_Name)
class = read.csv("ad_response/antidepressant_names_class.csv", header=T)[c(2,3)]
class = class[!duplicated(class),]

class_first_trial = merge(first_trial, class, by.x="Short_Name", by.y="Generic_name")

class_first_trial$drug_response = as.character(class_first_trial$drug_response)
class_first_trial$drug_response <- ifelse(class_first_trial$drug_response == "intermediate_multiple_antidepressant_and_antipsychotic" | class_first_trial$drug_response ==  "intermediate_antidepressant1_and_antipsychotic", "intermediate_with_antipsychotic", class_first_trial$drug_response)
class_first_trial$drug_response <- ifelse(class_first_trial$drug_response == "nonresponder_no_antipsychotics" | class_first_trial$drug_response == "nonresponder_with_antipsychotics", "nonresponder", class_first_trial$drug_response)

first_response_by_class_stats <- data.frame(table(class_first_trial$drug_claass, class_first_trial$drug_response))
write.csv(first_response_by_class_stats, paste0("first_drug_trial_response_by_class_stats_sample_sizes_",date,".csv"), row.names=T, quote=F)

## for demographics
### 1. sex
### 2. age at first ad
### 3. race
### 4. depression dx

# try new demo file
demo = fread("/data/davis_lab/shared/phenotype_data/biovu/delivered_data/sd_wide_pull/20230607_pull/20230607_sd_pull_person_all.txt.gz")
demo$DOB = as.Date(demo$birth_datetime, "%Y-%m-%d")
demo = demo[,c("GRID","gender_source_value","DOB" ,"race_source_value")]
dat = merge(class_first_trial, demo, by="GRID")
dat$first_drug_date = as.Date(dat$first_drug_date, "%Y-%m-%d")
dat$age_at_first_ad = as.numeric(dat$first_drug_date - dat$DOB)/365.25

# median age at first ad 
summary(dat$age_at_first_ad)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   # 0.00   31.88   47.96   46.70   61.16   89.85       1 


phecode = readRDS("/data/davis_lab/shared/phenotype_data/biovu/phecode_tables/medical_home_phecode_table_20210806_pull_remake_111822_2_distinct_dates_no_exclusions.Rds")
# 296.22 mdd
# 296.2 dep 
# 304 - adjustment rxn 
# 300.4 - dysthymic 

dep = phecode[c(1,487,488,512,499)]
colnames(dep) <- c("GRID","dep","mdd","adj_rxn","dysthymic")
dep_cases = subset(dep, dep=="TRUE" | mdd=="TRUE" | adj_rxn=="TRUE" | dysthymic=="TRUE")

get_demo = function(dat=dat, name=name){
	n_total = nrow(dat)
	n_female = nrow(subset(dat, gender_source_value=="F"))
	percent_female = round((n_female/nrow(dat))*100, digits=1)
	n_white = nrow(subset(dat, race_source_value=="W"))
	percent_white = round((n_white/nrow(dat))*100, digits=1)
	age = round(median(dat$age_at_first_ad, na.rm=T), digits=1)
	age1 = round(as.numeric(summary(dat$age_at_first_ad)[2]), digits=1)
	age3 = round(as.numeric(summary(dat$age_at_first_ad)[5]), digits=1)
	n_dep_cases = nrow(dat[(dat$GRID %in% dep_cases$GRID),])
	percent_dep_cases = round((n_dep_cases/nrow(dat))*100, digits=1)

	n_ssri = nrow(subset(dat, drug_claass=="SSRI"))
	n_snri = nrow(subset(dat, drug_claass=="SNRI"))
	n_tca = nrow(subset(dat, drug_claass=="TCA"))
	n_atypical = nrow(subset(dat, drug_claass=="Atypical"))
	n_maoi = nrow(subset(dat, drug_claass=="MAOI"))
	percent_ssri = round((n_ssri/nrow(dat))*100, digits=1)
	percent_snri = round((n_snri/nrow(dat))*100, digits=1)
	percent_tca = round((n_tca/nrow(dat))*100, digits=1)
	percent_atypical = round((n_atypical/nrow(dat))*100, digits=1)
	percent_maoi = round((n_maoi/nrow(dat))*100, digits=1)

	out = data.frame(name, n_total, n_female, percent_female, n_white, percent_white, age, age1, age3, n_dep_cases, percent_dep_cases, n_ssri, n_snri, n_tca, n_atypical, n_maoi, percent_ssri, percent_snri, percent_tca, percent_atypical, percent_maoi)
	return(out)
}

orig = dat

responder = get_demo(dat=subset(orig, drug_response=="responder"), name="responder")
nonresponder = get_demo(subset(orig, drug_response=="nonresponder"), name="nonresponder")
intermediate_antidepressant_and_antipsychotic = get_demo(subset(orig, drug_response=="intermediate_with_antipsychotic"), name="intermediate_with_antipsychotic")
intermediate_antidepressants_only = get_demo(subset(orig, drug_response=="intermediate_antidepressants_only"),  name="intermediate_antidepressants_only")

out = rbind(responder, intermediate_antidepressants_only, intermediate_antidepressant_and_antipsychotic, nonresponder)

out$n_female_percent = paste0(out$n_female, " (",out$percent_female, ")")
out$n_white_percent = paste0(out$n_white, " (",out$percent_white, ")")
out$median_age_iqr = paste0(out$age, " (",out$age1, " - ", out$age3, ")")
out$n_dep_percent = paste0(out$n_dep_cases, " (",out$percent_dep_cases, ")")

out$n_ssri_percent = paste0(out$n_ssri, " (",out$percent_ssri, ")")
out$n_snri_percent = paste0(out$n_snri, " (",out$percent_snri, ")")
out$n_tca_percent = paste0(out$n_tca, " (",out$percent_tca, ")")
out$n_atypical_percent = paste0(out$n_atypical, " (",out$percent_atypical, ")")
out$n_maoi_percent = paste0(out$n_maoi, " (",out$percent_maoi, ")")


write.csv(out, "first_trial_response_demographics_table.csv", quote=F, row.names=F)

## get time between first and last 
## find first drug trial response 
library(data.table)
setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/ad_response")
dat = fread("antidepressant_response_output_with_antipsychotics_052623_sd_pull_062123.txt", header=T)

first_date <- aggregate(MED_DATE ~ GRID, dat=dat, FUN=min)
last_date <- aggregate(MED_DATE ~ GRID, dat=dat, FUN=max)

dates = merge(first_date, last_date, by="GRID")

dates$months = interval(dates$MED_DATE.x, dates$MED_DATE.y) %/% days(1) / (365/12)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1.381   3.912  10.389  34.089  42.970 331.003 
