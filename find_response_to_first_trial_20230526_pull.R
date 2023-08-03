## find first drug trial response 
library(data.table)
setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/ad_response")
dat = fread("antidepressant_response_output_with_antipsychotics_052623_sd_pull_062123.txt", header=T)

first_episode <- aggregate(prescription_episode ~ GRID, dat=dat, FUN=min)
first_episode_dat <- merge(dat, first_episode, by=c("GRID", "prescription_episode"))
first_drug_first_episode <- setDT(first_episode_dat)[order(MED_DATE), head(.SD,1L),by="GRID"]
first_drug_first_episode2 <- first_drug_first_episode[, c("GRID","Short_Name", "MED_DATE")]
colnames(first_drug_first_episode2)[3] <- "first_drug_date"

first_drug_throughout_episode <- merge(first_drug_first_episode2, first_episode_dat, by=c("GRID","Short_Name"))
first_drug_final_response <- setDT(first_drug_throughout_episode)[order(-MED_DATE), head(.SD,1L),by="GRID"]
first_drug_final_response2 <- first_drug_final_response[, c("GRID","Short_Name","drug_response","first_drug_date","MED_DATE")]
colnames(first_drug_final_response2)[5] <- "last_drug_date"


# first_trial_out <- first_trial2(dat=output5_df_order)
date <- Sys.Date()
write.table(first_drug_final_response2, paste0("first_drug_trial_final_response_sd_pull_052623_",date,".txt"), row.names=F, quote=F, col.names=T, sep="\t")
