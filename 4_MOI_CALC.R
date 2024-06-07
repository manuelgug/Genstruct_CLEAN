
library(haven)
library(progress)
library(fs)
library(moire)

# Input: combined_df_merged_SAMPLING_only.RDS

# Outputs: all_samples_complete_filtered_MOIRE-RESULTS_SAMPLING_only_FOR_MOI.RDS (MOI); [province/region]_SAMPLING_MOIRE-RESULTS_FOR_ALLELE_FREQS.RDS (allele freqs and He for each population)

SAMPLING <- 2021 # 2021 or 2022

#######################################################
# calculate MOI and eMOI per run
#######################################################

combined_df_merged <- readRDS(paste0("combined_df_merged_", SAMPLING, "_only.RDS")) 

colnames(combined_df_merged)[1]<- c("sample_id")

# set MOIRE parameters
dat_filter <- moire::load_long_form_data(combined_df_merged)
burnin <- 1e4
num_samples <- 1e4
pt_chains <- seq(1, .5, length.out = 20)
pt_num_threads <- 20

# run moire
mcmc_results <- moire::run_mcmc(
  dat_filter, is_missing = dat_filter$is_missing,
  verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
  pt_chains = pt_chains, pt_num_threads = length(pt_chains),
  thin = 10); saveRDS(mcmc_results, paste0("all_samples_complete_filtered_MOIRE-RESULTS_", SAMPLING, "_only_FOR_MOI.RDS"))


#######################################################
# 8.- calculate He for provinces and regions
#######################################################

combined_df_merged <- readRDS(paste0("combined_df_merged_", SAMPLING, "_only.RDS"))

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
  saveRDS(mcmc_results, paste0(output_name, "_MOIRE-RESULTS_FOR_ALLELE_FREQS.RDS"))
}

# Create a list of data frames and corresponding years
data_frames <- list(combined_df_merged)
years <- list("2022")

# Loop over each province
for (i in seq_along(data_frames)) {
  year <- years[[i]]
  df <- data_frames[[i]]
  
  for (province in unique(df$province)) {
    province_df <- df[df$province == province, ]
    
    # Run MOIRE
    it_pr <-  paste0(province, "_", year)
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
    province_df <- df[df$region == region, ]
    
    # Run MOIRE
    it_re <-  paste0(region, "_", year)
    print(it_re)
    
    run_moire(province_df, it_re)
  }
}
