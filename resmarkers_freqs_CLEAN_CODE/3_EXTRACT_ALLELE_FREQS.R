library(haven)
library(progress)
library(fs)
library(moire)
library(dplyr)

SAMPLING <- 2022 # 2021 or 2022

#metadata
db <- read_dta('../DrugRes_2021_2022_DB_ALLDATA_1Jan2024.dta')

FINAL_MERGED_DIVERSITY_RESISTANCE_DF <- readRDS(paste0("FINAL_MERGED_DIVERSITY_RESISTANCE_DF", SAMPLING, ".RDS"))

#######################################################
# calculate MOI and eMOI per run
#######################################################

# 1) calculate heterozygosity of the population (He); pop = province, region
#import everything into lists
rds_files <- list.files(path = ".", pattern = paste0("\\MOIRE-RESULTS_FOR_ALLELE_FREQS_", SAMPLING,".RDS$"), full.names = TRUE)
rds_files <- rds_files[!rds_files %in% paste0("all_samples_complete_filtered_MOIRE-RESULTS_", SAMPLING, "only_FOR_MOI.RDS")] 


# Load each RDS file into the list with the file name as the list name
He_results_list <- list()

for (file in rds_files) {
  print(file)
  file_name <- tools::file_path_sans_ext(basename(file))
  He_results_list[[file_name]] <- readRDS(file)
}

# Loop through each element in He_results_list
processed_He_results <- data.frame()

for (i in seq_along(He_results_list)) {
  
  # Summarize He
  He_results <- moire::summarize_he(He_results_list[[i]])
  He_results$population <- names(He_results_list[i])
  
  # Add the processed He results to the list
  processed_He_results <- rbind(processed_He_results, He_results)
}

#formatting categories
processed_He_results$population <- gsub(paste0("_MOIRE-RESULTS_FOR_ALLELE_FREQS_",SAMPLING), "", processed_He_results$population)

# ######
# # keep amplicons with high He:
# 
# he_amps <- processed_He_results %>%
#   group_by(locus) %>%
#   summarize(mean = mean(post_stat_mean)) %>%
#   arrange(desc(mean))
# 
# #keep top n% amplicons with highest He (or use a numeric threshold like > 0.1 of He or whatever... test.)
# perc_25<- round(length(unique(he_amps$locus))*0.95)
# he_amps_top50 <- he_amps[1:perc_25,]
# 
# # FILTER
# processed_He_results <- processed_He_results[processed_He_results$locus %in% he_amps_top50$locus,]
# #####


library(stringr)
processed_He_results <- processed_He_results %>%
  mutate(geo = ifelse(str_detect(population, "North|South|Centre"), "region", "province"))

#processed_He_results$population <- gsub("TEST_|_202.*", "", processed_He_results$population)
#processed_He_results$year <- as.numeric(processed_He_results$year)

write.csv(processed_He_results, paste0("He_DIV_RES_", SAMPLING, ".csv"), row.names = F)




# Load each RDS file into the list with the file name as the list name
ALLELE_fREQ_results_list <- list()

for (file in rds_files) {
  print(file)
  file_name <- tools::file_path_sans_ext(basename(file))
  ALLELE_fREQ_results_list[[file_name]] <- readRDS(file)
}

# Loop through each element in ALLELE_fREQ_results_list
processed_ALLELE_fREQ_results <- data.frame()

for (i in seq_along(ALLELE_fREQ_results_list)) {
  
  # Summarize He
  ALLELE_fREQ_results <- moire::summarize_allele_freqs(ALLELE_fREQ_results_list[[i]])
  ALLELE_fREQ_results$population <- names(ALLELE_fREQ_results_list[i])
  
  # Add the processed He results to the list
  processed_ALLELE_fREQ_results <- rbind(processed_ALLELE_fREQ_results, ALLELE_fREQ_results)
}

#formatting categories
processed_ALLELE_fREQ_results$population <- gsub(paste0("_MOIRE-RESULTS_FOR_ALLELE_FREQS_", SAMPLING), "", processed_ALLELE_fREQ_results$population)

processed_ALLELE_fREQ_results <- processed_ALLELE_fREQ_results %>%
  mutate(geo = ifelse(str_detect(population, "North|South|Centre"), "region", "province"))

#processed_ALLELE_fREQ_results$population <- gsub("TEST_|_202.*", "", processed_ALLELE_fREQ_results$population)
#processed_ALLELE_fREQ_results$year <- as.numeric(processed_ALLELE_fREQ_results$year)

write.csv(processed_ALLELE_fREQ_results, paste0("ALLELE_fREQs_DIV_RES_", SAMPLING, ".csv"), row.names = F)

##########################################################################################################3
