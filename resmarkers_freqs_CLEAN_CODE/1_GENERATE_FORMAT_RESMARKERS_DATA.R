
SAMPLING <- 2022 # 2021 or 2022

path <- c("../")

combined_df_merged <- readRDS(paste0(path,"combined_df_merged_", SAMPLING, "_only.RDS")) 
combined_df_merged <- combined_df_merged[!(combined_df_merged$province %in% c("Maputo_Dry", "Manica_Dry")), ] # remove dry


db <- read_dta('../DrugRes_2021_2022_DB_ALLDATA_1Jan2024.dta')
colnames(db)
unique(db$run_id) #incomplete... must complete
unique(db$seasonality)


######################################################################
#-----------------------INPUT AND PREP DATA--------------------------#
######################################################################

# FILTERING RESULTS
SS <- length(unique(combined_df_merged$NIDA2))
cat("Final sample size is", as.character(SS))

LC <- length(unique(combined_df_merged$locus))
cat("Final Diversity loci count is", as.character(LC))

#allele formatting
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)
colnames(combined_df_merged)[1]<- c("SampleID")

#######################################################
#IMPORT MICROHAPS
#######################################################
# 1a.- complete run_id column with data from the actual runs

# Initialize an empty data frame to store the results
microhaps_data <- data.frame(NIDA2 = character(), NEW_run_id = character(), stringsAsFactors = FALSE)

directory_path <- "../results_v0.1.8_RESMARKERS_FIX/"

# Create a progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = length(dir_ls(path = directory_path, regexp = "_RESULTS_v0.1.8_FILTERED$"))
)

# Iterate through the folders ending with _RESULTS_v0.1.8
for (folder_path in dir_ls(path = directory_path, regexp = "_RESULTS_v0.1.8_FILTERED$")) {
  pb$tick()  # Update progress bar
  
  folder_name <- path_file(folder_path)
  file_path <- file.path(folder_path, "resmarker_microhap_table_global_max_0_filtered.csv")
  
  cat("\n")
  print(folder_name)
  
  # Read the contents of resmarker_microhap_table_global_max_0_filtered.csv
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
      microhaps_data <- rbind(microhaps_data, temp_df)
    }
  }
}

microhaps_data$NEW_run_id <- sub("_RESULTS_v0.1.8_FILTERED$", "", microhaps_data$NEW_run_id)

microhaps_data <- distinct(microhaps_data)

#this should be TRUE if the same number of nidas is in the db and the results from the grep
length(microhaps_data$NIDA2) == length(db$NIDA2) #it's not....
length(microhaps_data$NIDA2) - length(db$NIDA2) # 113 repeated nidas...

#check for repeated nidas 
repeated_nidas <- names(table(microhaps_data$NIDA2)[table(microhaps_data$NIDA2) > 1]) #there are duplicate nidas OOF
repeated_nidas_df<- microhaps_data[microhaps_data$NIDA2 %in% repeated_nidas,]
repeated_nidas_df <- repeated_nidas_df[order(repeated_nidas_df$NIDA2), ]
length(repeated_nidas_df$NIDA2)
length(unique(repeated_nidas_df$NIDA2))

#ask team about these nidas
#write.csv(repeated_nidas_df, paste0("repeated_nidas_",SAMPLING, ".csv"), row.names = F)


#######################################################
# 2.- genomic data from all runs (allele data, resmarkers, haplos?)
#######################################################
runs <- unique(paste0(microhaps_data$NEW_run_id, "_RESULTS_v0.1.8_FILTERED"))

folder_path <- paste0("../results_v0.1.8_RESMARKERS_FIX/", runs)
microhaps_files <- file.path(folder_path, "resmarker_microhap_table_global_max_0_filtered.csv")

# Import filtered allele data tables to the list
microhaps_list <- list()
for (file in microhaps_files) {
  microhaps <- read.csv(file)
  microhaps_list <- append(microhaps_list, list(microhaps))
}

#format the imported dfs
for (i in seq_along(microhaps_list)) {
  df <- microhaps_list[[i]]
  
  df$SampleID <- gsub("_S.*$", "", df$SampleID)
  df$SampleID <- gsub("_", ".", df$SampleID)
  df$SampleID <- gsub("-", ".", df$SampleID)
  df$SampleID <- gsub("N", "", df$SampleID)
  df$SampleID <- gsub("\\.0", "", df$SampleID)
  
  colnames(df)[1] <- "SampleID"
  
  df <- df %>% ### TEHERE IS A BUG WITH THE MASK THAT GENERATES MORE ALLELES THAN THERE REALLY ARE. THIS SNIPPET OF CODE COLLAPSES REPETITIONSAND SUMS THE READS AND FREQS. CRITICAL!!!
    group_by(SampleID, GeneID, Gene, MicrohapIndex, Microhaplotype) %>%
    summarize(reads = sum(Reads),
              norm.reads.locus = sum(norm.reads.locus)) %>%
    mutate(allele = paste(MicrohapIndex, ".", row_number(), sep = ""))
  
  
  # Update the modified data frame back in the list
  df<- as.data.frame(df)
  microhaps_list[[i]] <- df
}

names(microhaps_list) <- runs

#get rid of replicate nidas keep the sample with the most reads across runs for each one of the replicates
sum_reads <- function(df) {
  aggregate(reads ~ SampleID, data = df, sum)
}

summed_reads <- lapply(microhaps_list, sum_reads)

# Change colnames of reads for the name of the respective df
summed_reads <- lapply(names(summed_reads), function(df_name) {
  df <- summed_reads[[df_name]]
  names(df)[2] <- df_name
  return(df)
})

# Merge all data frames by sample_id
merged_df_dups <- Reduce(function(x, y) merge(x, y, by = "SampleID", all = TRUE), summed_reads)

#keep rows with replicates
merged_df_dups <- merged_df_dups[rowSums(is.na(merged_df_dups)) < length(colnames(merged_df_dups)) - 2, ] #keep rows with more than 14 NAs (that is, that have reads in more than one run)

# Exclude the first column and find the column with the maximum value for each row
merged_df_dups$BEST_RUN <- colnames(merged_df_dups)[apply(merged_df_dups[-1], 1, which.max) + 1]

# #did i found all replicates?
# if (all(repeated_nidas_df$NIDA2 %in% merged_df_dups$sample_id)) {
#   print("You found all replicates. Proceed with removal")
# }else{
#   "grab a coffee"
# }

#remove replicates, keep the best
for (i in 1:nrow(merged_df_dups)) {
  best_run <- merged_df_dups$BEST_RUN[i]
  SampleID <- merged_df_dups$SampleID[i]
  
  # Loop through microhaps_list
  for (j in seq_along(microhaps_list)) {
    df <- microhaps_list[[j]]
    
    # Exclude the df named after BEST_RUN
    if (names(microhaps_list)[j] != best_run) {
      microhaps_list[[j]] <- df[!(df$SampleID %in% SampleID), ]
    }
  }
}

#final check:
summed_reads <- lapply(microhaps_list, sum_reads)
summed_reads <- lapply(names(summed_reads), function(df_name) {
  df <- summed_reads[[df_name]]
  names(df)[2] <- df_name
  return(df)
})
merged_df_dups <- Reduce(function(x, y) merge(x, y, by = "SampleID", all = TRUE), summed_reads)
merged_df_dups <- merged_df_dups[rowSums(is.na(merged_df_dups)) < length(colnames(merged_df_dups)) - 2, ] #keep rows with more than 14 NAs (that is, that have reads in more than one run)

if (dim(merged_df_dups)[1] == 0){
  print("NO MORE REPLICATES.")
}else{
  "grab a coffee"
}

# Define a function to filter rows based on criteria: remove controls, basically
filter_rows <- function(df) {
  # Filter rows based on criteria
  filtered_df <- df[!grepl("^[A-Za-z]|3D", df$SampleID), ]
  return(filtered_df)
}

# Apply the filter_rows function to each dataframe in microhaps_list
microhaps_list <- lapply(microhaps_list, filter_rows)

# concat all dataframes together
combined_df_microhaps <- bind_rows(microhaps_list)

# calculate n.alleles for each locus of each sample if not done already during contaminant filtering
if (!("n.alleles" %in% colnames(combined_df_microhaps))){
  combined_df_microhaps <- combined_df_microhaps %>%
    group_by(SampleID, GeneID, MicrohapIndex) %>%
    mutate(n.alleles = n_distinct(allele))
}

# merge with metadata
colnames(combined_df_microhaps)[1]<- c("NIDA2")
combined_df_merged_microhaps <- merge(combined_df_microhaps, db[c("NIDA2", "year", "province", "region", "seasonality", "run_id")], by="NIDA2", all.y =T) #forcing adding all db nidas


################################# KEEP NIDA FROM CHOSEN SAMPLING #############################

combined_df_merged_microhaps <- combined_df_merged_microhaps[combined_df_merged_microhaps$NIDA2 %in% combined_df_merged$SampleID, ]

#saveRDS(combined_df_merged_microhaps, paste0("microhap_data_list", SAMPLING, ".RDS"))

##############################################################################################



#######################################################
# 3.- GENOMIC + DB MERGING, FILTERING ETC. (DATA PREP)
#######################################################

#allele_data_list <- readRDS(paste0("microhap_data_list", SAMPLING, ".RDS"))

# concat all dataframes together
combined_df <- bind_rows(allele_data_list)

# calculate n.alleles for each locus of each sample if not done already during contaminant filtering
if (!("n.alleles" %in% colnames(combined_df))){
  combined_df <- combined_df %>%
    group_by(sample_id, paste0(Gene, "_", MicrohapIndex)) %>%
    mutate(n.alleles = n_distinct(alleles))
}


if( sum(!(combined_df$NIDA2 %in% db$NIDA2)) == 0){
  print("All nidas in combined_merged_df are also the metadata db. No weird samples ✔")
}else{
  print("grab another coffee.")
}


#combined_df_merged <- combined_df


## FURTHER FILTERING
# 1) MAF filering (< 0.01)
combined_df <- combined_df[combined_df$norm.reads.locus  > 0.01, ]

# 3) remove bad diversity loci: those that are not in at least 100 samples with a read depth of 100. 
# Group by NIDA2 and locus, then summarize the total reads
locus_read_depth <- combined_df %>%
  group_by(NIDA2, locus = paste0(MicrohapIndex)) %>%
  summarize(read_depth = sum(reads))

# KEEP LOCI found in samples >=100 with read depth greater than 100
locus_counts <- locus_read_depth %>%
  group_by(locus) %>%
  summarize(samples_above_100 = sum(read_depth >= 100))

locus_counts <- na.omit(locus_counts)

loci_to_keep <- locus_counts[locus_counts$samples_above_100 >= 100,]$locus

combined_df <- combined_df[combined_df$MicrohapIndex %in% loci_to_keep, ]

library(forcats)

# MANAGE SEASONALITY
seasonality_factor <- as_factor(combined_df$seasonality)
seasonality_recode <- recode(seasonality_factor, "1" = "Rainy", "2" = "Dry")
seasonality_char <- as.character(seasonality_recode)
combined_df$seasonality <- seasonality_char

combined_df$province <- ifelse(
  is.na(combined_df$seasonality),
  combined_df$province,
  paste0(combined_df$province, "_", combined_df$seasonality)
)

# # JUST IN CASE... recalculate n.alleles for each locus of each sample
# combined_df <- combined_df %>%
#     group_by(NIDA2, locus) %>%
#     mutate(n.alleles = n_distinct(pseudo_cigar))

#recount n.alleles after filtering
combined_df <- combined_df %>%
  group_by(NIDA2, MicrohapIndex) %>%
  mutate(n.alleles = n_distinct(allele))

combined_df$province <- gsub(" ", "_", combined_df$province) # for cabo delgado

combined_df <- as.data.frame(combined_df)


# FILTERING RESULTS
SS <- length(unique(combined_df$NIDA2))
cat("Final sample size is", as.character(SS), " although some will be filtered out depending on the diversity results!!")

LC <- length(unique(combined_df$MicrohapIndex))
cat("Final Diversity loci count is", as.character(LC))

#save allele_data_list
#saveRDS(combined_df_microhaps, paste0("combined_df_microhaps_", SAMPLING, ".RDS"))




#######################################################
# Merge allele data and microhap data into a single df
#######################################################

#combined_df_merged_microhaps <- readRDS(paste0("combined_df_merged_microhaps_", SAMPLING, ".RDS"))

#rename microhap loci
combined_df_merged_microhaps$MicrohapIndex <- paste0(combined_df_merged_microhaps$Gene, "__", combined_df_merged_microhaps$MicrohapIndex)
#rename microhap alleles
combined_df_merged_microhaps$Microhaplotype <- paste0(combined_df_merged_microhaps$MicrohapIndex, "__", combined_df_merged_microhaps$Microhaplotype)

combined_df_merged_microhaps_SUBSET <- combined_df_merged_microhaps[c("NIDA2", "MicrohapIndex", "Microhaplotype", "reads", "province", "region", "run_id")]
combined_df_merged_SUBSET <- combined_df_merged[c("SampleID", "locus", "allele", "reads", "province", "region", "run_id")]

#rename columns
colnames(combined_df_merged_microhaps_SUBSET) <- c("sample_id", "locus", "allele", "reads", "province", "region", "run_id")
colnames(combined_df_merged_SUBSET) <- c("sample_id", "locus", "allele", "reads", "province", "region", "run_id")

#ids to extract from combined_df_merged_microhaps_SUBSET
ids_to_extract <- unique(combined_df_merged_SUBSET$sample_id)

#filter out sampleids not present in allele data
combined_df_merged_microhaps_SUBSET <- combined_df_merged_microhaps_SUBSET[combined_df_merged_microhaps_SUBSET$sample_id %in%ids_to_extract,]

#check sample_id correspondence between dfs
if (length(intersect(unique(combined_df_merged_microhaps_SUBSET$sample_id), unique(combined_df_merged_SUBSET$sample_id))) == length(unique(combined_df_merged_SUBSET$sample_id))){
  print("SampleIDs from allele and microhap data are the same")
}else{
  "grab a coffee."
}


FINAL_MERGED_DIVERSITY_RESISTANCE_DF <- rbind(combined_df_merged_SUBSET, combined_df_merged_microhaps_SUBSET)

FINAL_MERGED_DIVERSITY_RESISTANCE_DF$province <- gsub(" ", "_", FINAL_MERGED_DIVERSITY_RESISTANCE_DF$province) # for cabo delgado


if (SAMPLING == 2022){
  #need to remove "Maputo" and "Manica" alone because they are already with the *Dry and *Rainy suffixes in the data.. this happened due to the inclusion of the microhaplotype data. all good
  FINAL_MERGED_DIVERSITY_RESISTANCE_DF <- FINAL_MERGED_DIVERSITY_RESISTANCE_DF[FINAL_MERGED_DIVERSITY_RESISTANCE_DF$province != "Maputo",]
  FINAL_MERGED_DIVERSITY_RESISTANCE_DF <- FINAL_MERGED_DIVERSITY_RESISTANCE_DF[FINAL_MERGED_DIVERSITY_RESISTANCE_DF$province != "Manica",]
}

unique(FINAL_MERGED_DIVERSITY_RESISTANCE_DF$province)

sample_size_provinces <- FINAL_MERGED_DIVERSITY_RESISTANCE_DF %>%
  group_by(province) %>%
  summarise(unique_NIDA2_count = n_distinct(sample_id))

sample_size_provinces
sum(sample_size_provinces$unique_NIDA2_count)


combined_df_merged_nodry <- FINAL_MERGED_DIVERSITY_RESISTANCE_DF %>%
  filter(!grepl("Dry", province))

sample_size_regions <- combined_df_merged_nodry %>%
  group_by(region) %>%
  summarise(unique_NIDA2_count = n_distinct(sample_id))

sample_size_regions
sum(sample_size_regions$unique_NIDA2_count)


# # Replace contents of allele column based on condition
# FINAL_MERGED_DIVERSITY_RESISTANCE_DF <- FINAL_MERGED_DIVERSITY_RESISTANCE_DF %>%
#   mutate(allele = ifelse(grepl("^\\d", locus), paste0(locus, "___", allele), allele))


saveRDS(FINAL_MERGED_DIVERSITY_RESISTANCE_DF, paste0("FINAL_MERGED_DIVERSITY_RESISTANCE_DF", SAMPLING, ".RDS"))
