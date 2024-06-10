
library(haven)
library(progress)
library(fs)
library(moire)
library(dplyr)


SAMPLING <- 2022

FINAL_MERGED_DIVERSITY_RESISTANCE_DF <- readRDS(paste0("FINAL_MERGED_DIVERSITY_RESISTANCE_DF", SAMPLING, ".RDS"))


# #######################################################
# # calculate MOI and eMOI for each run
# #######################################################
# 
# # set MOIRE parameters
# dat_filter <- moire::load_long_form_data(FINAL_MERGED_DIVERSITY_RESISTANCE_DF)
# burnin <- 1e4
# num_samples <- 1e4
# pt_chains <- seq(1, .5, length.out = 20)
# pt_num_threads <- 20
# 
# # run moire
# mcmc_results <- moire::run_mcmc(
#   dat_filter, is_missing = dat_filter$is_missing,
#   verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
#   pt_chains = pt_chains, pt_num_threads = length(pt_chains),
#   thin = 10); saveRDS(mcmc_results, "all_samples_complete_filtered_MOIRE-RESULTS_", SAMPLING, "_only_FOR_MOI_DIV_RES.RDS")
# 


#######################################################
# calculate allele freqs for each population (per region/province)
#######################################################

# RUN IN CLUSTER:
# Define function to run MOIRE and save results
run_moire <- function(df, output_name) {
  colnames(df)[1] <- "sample_id"
  
  # set MOIRE parameters
  dat_filter <- moire::load_long_form_data(df)
  burnin <- 1e4
  num_samples <- 1e4
  pt_chains <- seq(1, .5, length.out = 20)
  
  # run moire
  mcmc_results <- moire::run_mcmc(
    dat_filter, is_missing = dat_filter$is_missing,
    verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
    pt_chains = pt_chains, pt_num_threads = length(pt_chains),
    thin = 10
  )
  
  # checkpoint
  saveRDS(mcmc_results, paste0(output_name, paste0("_MOIRE-RESULTS_FOR_ALLELE_FREQS_", SAMPLING, ".RDS")))
}

# Create a list of data frames and corresponding years
data_frames <- list(FINAL_MERGED_DIVERSITY_RESISTANCE_DF)
years <- list(SAMPLING)

# Loop over each province
for (i in seq_along(data_frames)) {
  year <- years[[i]]
  df <- data_frames[[i]]
  
  for (province in unique(df$province)) {
    province_df <- df[df$province == province, ]
    
    # Run MOIRE
    it_pr <-  paste0(province, "_RES_DIV_", year)
    print(it_pr)
    
    run_moire(province_df, it_pr)
  }
}

data_frames[[1]] <- data_frames[[1]][!(data_frames[[1]]$province %in% c("Maputo_Dry", "Manica_Dry")), ] # remove DRY season pops from region analysis

# Loop over each region
for (i in seq_along(data_frames)) {
  year <- years[[i]]
  df <- data_frames[[i]]
  
  for (region in unique(df$region)) {
    region_df <- df[df$region == region, ]
    
    # Run MOIRE
    it_re <-  paste0(region, "_RES_DIV_", year)
    print(it_re)
    
    run_moire(region_df, it_re)
  }
}