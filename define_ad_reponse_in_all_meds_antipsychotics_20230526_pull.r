setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230511/ad_response/")
# Rscript define_ad_response_in_all_meds_antipsychotics_20230607.r

library(data.table)
library(dplyr)
library(tidyverse)

meds_clean = fread("sd_wide_antidepressants_from_052623_pull.csv")
meds_clean$matched_drug = tolower(meds_clean$matched_drug)
meds_clean = subset(meds_clean, as.Date(meds_clean$drug_exposure_start_date, "%Y-%m-%d") > as.Date("1990-01-01"))

ap = fread("/data/davis_lab/allie/atap_weight/data/sd_wide_atypical_antipsychotics_dbc_from_052623_pull.csv", header=T)
ap2 <- ap[,c("GRID","matched_drug_generic","drug_exposure_start_date")]

ap2_u <- ap2[!(duplicated(ap2)),]
ap_in_ad = ap2_u[(ap2_u$GRID %in% meds_clean$GRID),]

# combine antidepressants and antipsychotics 
# meds_ad2 = epic_meds_ad[,c(15,16,17,3)]
meds_ad = meds_clean[,c(1,6,8)]
colnames(meds_ad) <- c("GRID","Short_Name","MED_DATE")
meds_ad$drug_type = "antidepressant"
# colnames(meds_ad2)[3] <- "drug_class"
meds_ad_u = meds_ad[!duplicated(meds_ad),]

colnames(ap_in_ad) <- c("GRID","Short_Name","MED_DATE")
# ap_in_ad$drug_class <- "antipsychotic"
ap_in_ad$drug_type <- "antipsychotic"
ap_in_ad$Short_Name = tolower(ap_in_ad$Short_Name)

ad_ap <- rbind(meds_ad_u, ap_in_ad)

### rerun ad response using epic antidepressants + antipsychotics 

ad_ap$MED_DATE <- as.Date(ad_ap$MED_DATE, "%Y-%m-%d")
ad_ap_order <- ad_ap[order(ad_ap$GRID, ad_ap$MED_DATE),]
rownames(ad_ap_order) <- seq(1:nrow(ad_ap_order))

## build rx episodes from antidepressants
## apply rx episodes to antipsychotics
## require ad to have >6 weeks
## require ap to have >2 weeks

find_rx_epidodes = function(input=input){
    output1 <- input %>%
        group_by(GRID) %>%
        arrange(MED_DATE) %>%
        mutate(prev_drug_date = dplyr::lag(MED_DATE, n = 1, default = NA)) %>%
        mutate(diff_weeks_drug = as.numeric(difftime(MED_DATE, prev_drug_date, units = "weeks"))) %>%
        mutate(prescription_episode = ifelse(diff_weeks_drug > 26, seq_along(diff_weeks_drug), 1)) %>%
        mutate(prescription_episode = ifelse(is.na(prescription_episode), 999, prescription_episode)) %>%
        mutate(prescription_episode = ifelse(prescription_episode == 1, NA, prescription_episode)) %>%
        fill(prescription_episode) %>%
        mutate(prescription_episode = ifelse(prescription_episode == 999, 1, prescription_episode))
    return(output1)
    }

# require each med to have evidence of being prescribed for at least 6 weeks
require_adequate_trial <- function(input=input, min_weeks=min_weeks){
    output2 <- input %>%
        group_by(GRID, prescription_episode, Short_Name) %>%
        arrange(MED_DATE) %>%
        mutate(prev_same_drug_date = dplyr::lag(MED_DATE, n = 1, default = NA)) %>%
        mutate(diff_weeks_same_drug = as.numeric(difftime(MED_DATE, prev_same_drug_date, units = "weeks"))) %>%
        mutate(prescription_episode_same_drug = sum(as.numeric(diff_weeks_same_drug),na.rm=T)) %>%
        mutate(adequate_prescription_episode_same_drug = ifelse(prescription_episode_same_drug >= min_weeks, 1, 0))
    output3 <- subset(output2, adequate_prescription_episode_same_drug=="1")
    return(output3)
    }

apply_rx_episodes <- function(ad_rx_episodes = ad_rx_episodes, ap_input=ap_input){
    rx_ep_max = aggregate(MED_DATE ~ GRID+prescription_episode, data=ad_rx_episodes, FUN=max)
    rx_ep_min = aggregate(MED_DATE ~ GRID+prescription_episode, data=ad_rx_episodes, FUN=min)
    colnames(rx_ep_max)[3] <- 'max_med_date'
    colnames(rx_ep_min)[3] <- 'min_med_date'
    rx_ep_dates = merge(rx_ep_min, rx_ep_max, by=c("GRID", "prescription_episode"))
    rx_ep_dates_ap = merge(rx_ep_dates, ap_input, by="GRID")
    ap_with_rx_dates = subset(rx_ep_dates_ap, MED_DATE >= min_med_date & MED_DATE <= max_med_date)
    return(ap_with_rx_dates)
    }

## find switches
## find response 
## find ap adjunctive if responder or intermediate 

add_ap_flag <- function(ad_dat = dat_dat, ap_dat = ap_dat){
    dat = rbind(ad_adequate_trial, ap_adequate_trial)
    dat2 = dat %>%
        group_by(GRID, prescription_episode) %>%
        arrange(MED_DATE) %>%
        mutate(antipsychotic_present = case_when('antipsychotic' %in% drug_type ~ TRUE))
    out = subset(dat2, drug_type == "antidepressant")
    return(out)
    }

# drug switch = yes if the drug is not the last date in the episode and there are no later dates for the same drug
# drug switch = no if there are later dates for the same drug within the episode or if the date is the last in the episode 

# non responder = drug switch
    # 1. non-responder no antipsychotics = drug_switch == "yes" & is.na(antipsychotic_present)
    # 2. non-responder with antipsychotics = drug_switch == "yes" & antipsychotic_present == "TRUE"
# intermeidate = 
    # 1. no drug switch + multiple antidepressants within the episode 
    # 2. no drug switch + antipsychotic within the episode 
    # 3. no drug switch + 1 drug in episode + antipsychotic flag 
# responder = no drug switch + 1 drug in episode + no antipsychotic flag

# intermediate_antidepressants_only = drug_switch == "no" & drug_distinct > 1 & is.na(antipsychotic_present)
# intermediate_antidepressant1_and_antipsychotic = drug_switch == "no" & drug_distinct == 1 & antipsychotic_present == "TRUE"
# intermediate_multiple_antidepressants_and_antipsychotic = drug_switch == "no" & drug_distinct > 1 & antipsychotic_present == "TRUE"

find_response_original <- function(input=input){
    # find next date of same drug to determine switch
    output4 <- input %>%
        group_by(GRID, prescription_episode, Short_Name) %>%
        arrange(MED_DATE) %>%
        mutate(next_same_drug_date = dplyr::lead(MED_DATE, n = 1, default = NA))
    # define switch and drug response 
    output5 <- output4 %>%
        group_by(GRID, prescription_episode) %>%
        arrange(MED_DATE) %>%
        mutate(drug_distinct = n_distinct(Short_Name,na.rm=T)) %>% ## change this to be n_drugs on (or after?) the same date
        mutate(last_drug_in_episode_prescr_date = max(as.Date(MED_DATE, "%Y-%m-%d"))) %>%
        mutate(last_drug_in_episode_prescr = ifelse(as.Date(MED_DATE, "%Y-%m-%d") == last_drug_in_episode_prescr_date, 'yes', 'no')) %>%
        mutate(drug_switch = ifelse(!is.na(next_same_drug_date), "no", 
            ifelse(is.na(next_same_drug_date) & last_drug_in_episode_prescr == "yes", "no", "yes"))) %>%
        mutate(drug_response = ifelse(drug_switch == "no" & drug_distinct == "1" & is.na(antipsychotic_present), "responder",
        ifelse(drug_switch == "yes" & is.na(antipsychotic_present), "nonresponder_no_antipsychotics", 
            ifelse(drug_switch == "yes" & antipsychotic_present == "TRUE", "nonresponder_with_antipsychotics",
                ifelse(drug_switch == "no" & drug_distinct > 1 & is.na(antipsychotic_present), "intermediate_antidepressants_only",
                    ifelse(drug_switch == "no" & drug_distinct == "1" & antipsychotic_present == "TRUE", "intermediate_antidepressant1_and_antipsychotic",
                        ifelse(drug_switch == "no" & drug_distinct > 1 & antipsychotic_present == "TRUE", "intermediate_multiple_antidepressant_and_antipsychotic", "NA")))))))
    output5_df <- as.data.frame(output5)
    output5_df_order <- output5_df[order(output5_df$GRID,output5_df$MED_DATE),]
    return(output5_df_order)
    }



## run all functions
ad_input = subset(ad_ap_order, drug_type=="antidepressant")
ad_rx_episodes = find_rx_epidodes(input=ad_input)
ad_adequate_trial = require_adequate_trial(input=ad_rx_episodes, min_weeks = 6)

ap_input = subset(ad_ap_order, drug_type=="antipsychotic")
ap_rx_episodes = apply_rx_episodes(ad_rx_episodes = ad_rx_episodes, ap_input=ap_input)
ap_rx_episodes = ap_rx_episodes[c(1,2,5:7)]
ap_adequate_trial = require_adequate_trial(input=ap_rx_episodes, min_weeks=2)

ad_with_ap_flag = add_ap_flag(ad_dat = ad_adequate_trial, ap_dat = ap_adequate_trial)

response = find_response_original(input=ad_with_ap_flag)
# response1 = find_response1(input=ad_with_ap_flag)
# response2 = find_response2(input=ad_with_ap_flag)

print("n response rows")
print(nrow(response))
# [1] 7,152,656

print("n response samples")
print(nrow(response[!duplicated(response$GRID),]))
# [1] 315935

rownames(response) <- NULL
write.table(response, "antidepressant_response_output_with_antipsychotics_052623_sd_pull_062123.txt", col.names=T, sep="\t", row.names=F, quote=F)

print("response numbers")
print(summary(as.factor(response$drug_response)))
#                intermediate_antidepressant1_and_antipsychotic 
#                                                 449176 
#                      intermediate_antidepressants_only 
#                                                2419935 
# intermediate_multiple_antidepressant_and_antipsychotic 
#                                                 936862 
#                         nonresponder_no_antipsychotics 
#                                                  66463 
#                       nonresponder_with_antipsychotics 
#                                                  18544 
#                                              responder 
#                                                3261676



                                               