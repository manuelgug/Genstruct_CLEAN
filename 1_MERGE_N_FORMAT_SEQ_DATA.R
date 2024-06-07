

## MERGE AND FORMAT SEQ DATA. 

# Input: seq run folders after pipeline and filtering
# Output: allele_data_list.RDS


library(haven)
library(progress)
library(fs)
library(ggplot2)
library(dplyr)


#######################################################
# metadata (provided by Simone)
#######################################################

db <- read_dta('DrugRes_2021_2022_DB_ALLDATA_1Jan2024.dta')
colnames(db)
unique(db$run_id) #incomplete... must complete


#######################################################
# merge and format
#######################################################

# Initialize an empty data frame to store the results
result_df <- data.frame(NIDA2 = character(), NEW_run_id = character(), stringsAsFactors = FALSE)

directory_path <- "results_v0.1.8_RESMARKERS_FIX/"

# Create a progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = length(dir_ls(path = directory_path, regexp = "_RESULTS_v0.1.8$"))
)

# Iterate through the folders ending with _RESULTS_v0.1.8
for (folder_path in dir_ls(path = directory_path, regexp = "_RESULTS_v0.1.8$")) {
  pb$tick()  # Update progress bar
  
  folder_name <- path_file(folder_path)
  file_path <- file.path(folder_path, "amplicon_coverage.txt")
  
  cat("\n")
  print(folder_name)
  
  # Read the contents of /quality_report/amplicon_stats.txt
  if (file.exists(file_path)) {
    sample_coverage_content <- readLines(file_path)
    
    # Truncate each line after the first tab (\t) and return unique values and edit nida format to match the one from the db
    truncated_values <- unique(sapply(strsplit(sample_coverage_content, "\t"), function(x) x[1]))
    truncated_values <- gsub("_S.*$", "", truncated_values)
    truncated_values <- gsub("_", ".", truncated_values)
    truncated_values <- gsub("-", ".", truncated_values)
    truncated_values <- gsub("N", "", truncated_values)
    
    # Extract NIDA2 from runs using grep
    nida_values <- grep(paste(db$NIDA2, collapse = "|"), truncated_values, value = TRUE) #not all samples from runs are needed, only those in the db, hence this
    
    # Create a data frame with NIDA2 and folder_name
    if (length(nida_values) > 0) {
      temp_df <- data.frame(NIDA2 = nida_values, NEW_run_id = folder_name, stringsAsFactors = FALSE)
      
      # Append the results to the main data frame
      result_df <- rbind(result_df, temp_df)
    }
  }
}

result_df$NEW_run_id <- sub("_RESULTS_v0.1.8$", "", result_df$NEW_run_id)

#this should be TRUE if the same number of nidas is in the db and the results from the grep
length(result_df$NIDA2) == length(db$NIDA2) #it's not....
length(result_df$NIDA2) - length(db$NIDA2) # 113 repeated nidas...

#check for repeated nidas 
repeated_nidas <- names(table(result_df$NIDA2)[table(result_df$NIDA2) > 1]) #there are duplicate nidas OOF
repeated_nidas_df<- result_df[result_df$NIDA2 %in% repeated_nidas,]
repeated_nidas_df <- repeated_nidas_df[order(repeated_nidas_df$NIDA2), ]
length(repeated_nidas_df$NIDA2)
length(unique(repeated_nidas_df$NIDA2))

#ask team about these nidas
write.csv(repeated_nidas_df, "repeated_nidas.csv", row.names = F)


#######################################################
# genomic data from all runs 
#######################################################
runs <- unique(paste0(result_df$NEW_run_id, "_RESULTS_v0.1.8_FILTERED"))

folder_path <- paste0("results_v0.1.8_RESMARKERS_FIX/", runs)
allele_data_files <- file.path(folder_path, "allele_data_global_max_0_filtered.csv")

# Import filtered allele data tables to the list
allele_data_list <- list()
for (file in allele_data_files) {
  allele_data <- read.csv(file)
  allele_data_list <- append(allele_data_list, list(allele_data))
}

#format the imported dfs
for (i in seq_along(allele_data_list)) {
  df <- allele_data_list[[i]]
  
  df$sampleID <- gsub("_S.*$", "", df$sampleID)
  df$sampleID <- gsub("_", ".", df$sampleID)
  df$sampleID <- gsub("-", ".", df$sampleID)
  df$sampleID <- gsub("N", "", df$sampleID)
  df$sampleID <- gsub("\\.0", "", df$sampleID)
  
  colnames(df)[1] <- "sample_id"
  
  df <- df %>% ### TEHERE IS A BUG WITH THE MASK THAT GENERATES MORE ALLELES THAN THERE REALLY ARE. THIS SNIPPET OF CODE COLLAPSES REPETITIONSAND SUMS THE READS AND FREQS. CRITICAL!!!
    group_by(sample_id, locus, pseudo_cigar) %>%
    summarize(reads = sum(reads),
              norm.reads.locus = sum(norm.reads.locus),
              Category = first(Category)) %>%
    mutate(allele = paste(locus, ".", row_number(), sep = ""))
  
  
  # Update the modified data frame back in the list
  df<- as.data.frame(df)
  allele_data_list[[i]] <- df
}

names(allele_data_list) <- runs

#get rid of replicate nidas keep the sample with the most reads across runs for each one of the replicates
sum_reads <- function(df) {
  aggregate(reads ~ sample_id, data = df, sum)
}

summed_reads <- lapply(allele_data_list, sum_reads)

# Change colnames of reads for the name of the respective df
summed_reads <- lapply(names(summed_reads), function(df_name) {
  df <- summed_reads[[df_name]]
  names(df)[2] <- df_name
  return(df)
})

# Merge all data frames by sample_id
merged_df_dups <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), summed_reads)

#keep rows with replicates
merged_df_dups <- merged_df_dups[rowSums(is.na(merged_df_dups)) < length(colnames(merged_df_dups)) - 2, ] #keep rows with more than 14 NAs (that is, that have reads in more than one run)

# Exclude the first column and find the column with the maximum value for each row
merged_df_dups$BEST_RUN <- colnames(merged_df_dups)[apply(merged_df_dups[-1], 1, which.max) + 1]

#remove replicates, keep the best
for (i in 1:nrow(merged_df_dups)) {
  best_run <- merged_df_dups$BEST_RUN[i]
  sample_id <- merged_df_dups$sample_id[i]
  
  # Loop through allele_data_list
  for (j in seq_along(allele_data_list)) {
    df <- allele_data_list[[j]]
    
    # Exclude the df named after BEST_RUN
    if (names(allele_data_list)[j] != best_run) {
      allele_data_list[[j]] <- df[!(df$sample_id %in% sample_id), ]
    }
  }
}

#final check:
summed_reads <- lapply(allele_data_list, sum_reads)
summed_reads <- lapply(names(summed_reads), function(df_name) {
  df <- summed_reads[[df_name]]
  names(df)[2] <- df_name
  return(df)
})
merged_df_dups <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), summed_reads)
merged_df_dups <- merged_df_dups[rowSums(is.na(merged_df_dups)) < length(colnames(merged_df_dups)) - 2, ] #keep rows with more than 14 NAs (that is, that have reads in more than one run)

if (dim(merged_df_dups)[1] == 0){
  print("NO MORE REPLICATES.")
}else{
  "grab a coffee"
}

# Define a function to filter rows based on criteria: remove controls, basically
filter_rows <- function(df) {
  # Filter rows based on criteria
  filtered_df <- df[!grepl("^[A-Za-z]|3D", df$sample_id), ]
  return(filtered_df)
}

# Apply the filter_rows function to each dataframe in allele_data_list
allele_data_list <- lapply(allele_data_list, filter_rows)

#save allele_data_list
saveRDS(allele_data_list, "allele_data_list.RDS")