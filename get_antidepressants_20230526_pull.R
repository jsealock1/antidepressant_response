# setwd("/data/davis_lab/allie/atap_weight/")
setwd("/data/davis_lab/sealockj/projects/mdd_wbc_long/ad_response/add_antipsychotics/sd_pull_20230607/ad_response/")

library(data.table)
library(stringr)

stringify <- function(xx, sep, add_quote=T) {
  str <- as.character(xx)
  if (add_quote) {str <- paste0("'",xx,"'")}
  str <-  paste0(str,collapse = sep)
  str
}

## load list of atypical antipsychotics
# list_atap <- fread("data/my_curated_atypical_antipsychotics_list_with_sulpiride_for_dbc_query_040323.csv")
list_atap <- fread("antidepressant_names_class.csv")
names(list_atap) <- gsub(" ","_",names(list_atap))
# list_atap <- list_atap[Generic_name!="sulpiride"]

list_atap_vec <- sort(toupper(c(list_atap$Generic_name, list_atap[!is.na(Brand_name) & Brand_name!="",Brand_name])))
list_atap_regex <-  stringify(list_atap_vec,"|", add_quote=F)

## load in drug_exposure table (takes ~20-25 mins with 36 threads on vgi01)
d <- fread("/data/davis_lab/shared/phenotype_data/biovu/delivered_data/sd_wide_pull/20230607_pull/20230607_sd_pull_drug_exposure_all.txt.gz",quote="")

## regex search for the relevant drugs
## this step takes ~1 hour
d_atap <- d[grepl(list_atap_regex,toupper(drug_source_value))]
d_atap[,matched_drug := str_extract(toupper(drug_source_value), list_atap_regex)]

## convert brand names to generic
list_atap_brand <- list_atap[!is.na(Brand_name),.(toupper(Generic_name), toupper(Brand_name))] 
names(list_atap_brand) <- c("generic","brand")
nrow(d_atap) # 4891490
d_atap <- merge(d_atap, list_atap_brand, by.x="matched_drug", by.y="brand", all.x=T)
nrow(d_atap) # 4891490 - just checking that rows in dat didn't match to multiple rows in list_atap_brand, creating duplicate entries
d_atap[,matched_drug_generic := ifelse(is.na(generic),matched_drug,generic)]

stopifnot(nrow(d_atap[is.na(matched_drug_generic)])==0)

## map drug_type_concept_id
drug_type_map <- data.frame(
  drug_type_concept_id = c(38000175, 38000177, 38000178, 38000179, 38000180, 44787730, 2002843132),
  drug_type_concept_name = c("Prescription dispensed in pharmacy", "Prescription written","Medication list entry",
                   "Physician administered drug (identified as procedure)", "Inpatient administration",
                   "Patient Self-Reported Medication", "CPT Codes in Drug")
)
d_atap <- merge(d_atap, drug_type_map, by="drug_type_concept_id", all.x=T)

## check column missingness
mean(d_atap$provider_id==0) # 69% missing
mean(d_atap$stop_reason=="") # 96% missing
mean(is.na(d_atap$refills)) # 97% missing
mean(is.na(d_atap$route_concept_id)) # 31% missing
mean(d_atap$route_source_value=="") # 29% missing
mean(d_atap$dose_unit_source_value=="") # 57% missing
mean(is.na(d_atap$x_src)) # 100% 
mean(is.na(d_atap$lot_number)) # 100%
mean(is.na(d_atap$quantity)) # 75% missing
mean(is.na(d_atap$days_supply)) # 99% missing
mean(d_atap$x_drug_type_source_concept_id==0) # 100% missing

## create final cleaned table
d_atap_clean <- d_atap[,.(GRID,drug_exposure_id,drug_type_concept_id,drug_type_concept_name,matched_drug,matched_drug_generic,drug_concept_id,drug_exposure_start_date,drug_exposure_end_date,quantity,route_concept_id,provider_id,visit_occurrence_id,drug_source_value,drug_source_concept_id,route_source_value,dose_unit_source_value)]
setorder(d_atap_clean,GRID,drug_exposure_start_date)
d_atap_clean[,matched_drug_generic := str_to_upper(matched_drug_generic)]
d_atap_clean[,matched_drug := str_to_upper(matched_drug)]

fwrite(d_atap_clean, "sd_wide_antidepressants_from_20230607_pull.csv")
length(unique(d_atap_clean$GRID)) # 156372 individuals

## get drug_exposure_id for mapping to x_drug_exposure table
ids <- unique(d_atap_clean$drug_exposure_id)

## load X_DRUG_EXPOSURE table, which contains some additional fields
## this takes ~30 minutes to load
xd <- fread("/data/davis_lab/shared/phenotype_data/biovu/delivered_data/sd_wide_pull/20230607_pull/20230607_sd_pull_x_drug_exposure_all.txt.gz",quote="")
xd_subset <- xd[drug_exposure_id %in% ids]
length(unique(xd_subset$GRID)) # 155422 individuals

## check column missingness
mean(xd_subset$x_doc_stype=="") # 53% missing
mean(xd_subset$x_doc_type=="") # 0% missing
mean(xd_subset$x_dose=="") # 46% missing
mean(xd_subset$x_drug_form=="") # 38% missing
mean(xd_subset$x_strength=="") # 12% missing
mean(xd_subset$x_ndc=="") # 63% missing
mean(xd_subset$x_frequency=="") # 38% missing
mean(xd_subset$x_quantity_unit=="") # 87% missing
mean(xd_subset$x_rxnorm_source=="") # <1% missing

xd_subset_clean <- xd_subset[,.(GRID,drug_exposure_id,x_doc_type,x_dose,x_drug_form,x_strength,x_frequency,x_ndc,x_rxnorm_source)]

## merge drug_exposure and x_drug_exposure tables
merged <- merge(d_atap_clean, xd_subset_clean, by=c("drug_exposure_id","GRID"),all.x=T)
nrow(merged[is.na(x_doc_type)]) # 24,199 med records did not map to the X table (<1%)

## write out merged file
# fwrite(merged, "data/sd_wide_atypical_antipsychotics_merged_xtable_dbc_from_052623_pull.csv")
fwrite(merged, "sd_wide_antidepressants_merged_xtable_dbc_from_20230607_pull.csv")

