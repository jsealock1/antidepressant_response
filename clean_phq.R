## phq 8 = 122,041 entries
## phq 2 = 120,577 entries
## phq 9 = 63 entries

setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/phq_data")

dat = fread("/data/davis_lab/shared/phenotype_data/biovu/delivered_data/sd_wide_pull/20230607_pull/20230607_sd_pull_x_survey.txt.gz", header=T)

phq2 = subset(dat, form_name=="MHAV PROM PHQ-2" & question=="PHQ-2 Score")
phq8 = subset(dat, form_name=="MHAV PROM PHQ-8" & question=="PHQ-8 SCORE" | question=="PHQ-8 Severity Score")



clean_phq = function(phq=phq, max_value=max_value){
	phq = subset(phq, answer<=max_value)
	phq$date = as.Date(phq$answer_datetime, "%Y-%m-%d")
	phq = phq[,c("GRID","answer","form_name","date")]
	phq$answer = as.numeric(phq$answer)
	phq_out = aggregate(answer ~ GRID + date, FUN=mean, data=phq)
	return(phq_out)
}

phq2_clean = clean_phq(phq=phq2, max_value=6)
phq8_clean = clean_phq(phq=phq8, max_value=24)

write.table(phq2_clean, "phq2_clean_data_sd_wide_pull_20230607.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(phq8_clean, "phq8_clean_data_sd_wide_pull_20230607.txt", col.names=T, row.names=F, quote=F, sep="\t")