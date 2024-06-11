
library(haven)
library(progress)
library(fs)
library(dplyr)


# Input: allele_data_list.RDS

# Output: combined_df_merged_SAMPLING_only.RDS

######################################################################
#-----------------------INPUT AND PREP DATA--------------------------#
######################################################################

#metadata
db <- read_dta('DrugRes_2021_2022_DB_ALLDATA_1Jan2024.dta')

# indicate sampling
SAMPLING <- 2022 # 2022 or 2021


#######################################################
# GENOMIC + DB MERGING, FILTERING ETC. (DATA PREP)
#######################################################

allele_data_list <- readRDS("allele_data_list.RDS")

# concat all dataframes together
combined_df <- bind_rows(allele_data_list)

# calculate n.alleles for each locus of each sample if not done already during contaminant filtering
if (!("n.alleles" %in% colnames(combined_df))){
  combined_df <- combined_df %>%
    group_by(sample_id, locus) %>%
    mutate(n.alleles = n_distinct(pseudo_cigar))
}

# merge with metadata
colnames(combined_df)[1]<- c("NIDA2")
combined_df_merged <- merge(combined_df, db[c("NIDA2", "year", "province", "region", "seasonality", "run_id")], by="NIDA2", all.y =T) #forcing adding all db nidas

#check for columns with NAs. THESE SAMPLES WERE BAD QUALITY AND THUS FILTERED OUT DURING CONTAMINANTS FILTERING
removed_samples <- combined_df_merged[is.na(combined_df_merged$Category),]
removed_samples$NIDA2

# Remove rows with NIDA2 matching values in removed_samples
combined_df_merged <- combined_df_merged[!combined_df_merged$NIDA2 %in% removed_samples$NIDA2, ]


#sanity check
if (sum(is.na(combined_df_merged[, !colnames(combined_df_merged) %in% "seasonality"])) == 0) {
  print("No NAs ✔")
} else {
  print("grab a coffee.")
}

if( sum(!(combined_df_merged$NIDA2 %in% db$NIDA2)) == 0){
  print("All nidas in combined_merged_df are also the metadata db. No weird samples ✔")
}else{
  print("grab another coffee.")
}


#KEEP ONLY DIVERTSITY LOCI!
# are there samples with no Diversity loci sequenced?
unique_categories <- combined_df_merged %>%
  group_by(NIDA2) %>%
  summarize(unique_categories = toString(unique(Category)))

more_removed_sample <- unique_categories[!grepl("Diversity", unique_categories$unique_categories), ] # 1 sample didn't have diversity loci, probably due to a missing pool or something... no problem.
more_removed_sample

combined_df_merged <- combined_df_merged[combined_df_merged$Category == "Diversity",] 

#KEEP ONLY SAMPLING samples!
combined_df_merged <- combined_df_merged[combined_df_merged$year == SAMPLING,] 


## FURTHER FILTERING

# 1) MAF filering (< 0.01)
combined_df_merged <- combined_df_merged[combined_df_merged$norm.reads.locus  > 0.01, ]

# 2) check for coverage (>50 loci with >threshold reads)
# Define thresholds for read depth
thresholds <- c(25, 50, 100, 200)
count_list <- list()

# Loop over each threshold
for (threshold in thresholds) {
  # Calculate unique loci counts for the current threshold
  count <- combined_df_merged %>%
    group_by(NIDA2, locus) %>%
    summarize(total_reads = sum(reads)) %>%
    group_by(NIDA2) %>%
    filter(total_reads > threshold) %>%
    summarize(!!paste("unique_loci_at_min_", threshold, "_reads", sep = "") := n_distinct(locus))

  count_list[[paste("count_", threshold, sep = "")]] <- count
}

result_df <- Reduce(function(x, y) left_join(x, y, by = "NIDA2"), count_list)

#keep samples with >50 unique loci >= 100 reads
result_df <- na.omit(result_df[result_df$unique_loci_at_min_100_reads > 50,][c(1,4)])
samples_to_keep <- result_df$NIDA2 

combined_df_merged <- combined_df_merged[combined_df_merged$NIDA2 %in% samples_to_keep, ]


# 3) remove bad diversity loci: those that are not in at least 100 samples with a read depth of 100. 
# Group by NIDA2 and locus, then summarize the total reads
locus_read_depth <- combined_df_merged %>%
  group_by(NIDA2, locus) %>%
  summarize(read_depth = sum(reads))

# Count the number of samples with read depth greater than 100 for each locus
locus_counts <- locus_read_depth %>%
  group_by(locus) %>%
  summarize(samples_above_100 = sum(read_depth > 100))

loci_to_keep <- locus_counts$locus

combined_df_merged <- combined_df_merged[combined_df_merged$locus %in% loci_to_keep, ]

library(forcats)

# MANAGE SEASONALITY
seasonality_factor <- as_factor(combined_df_merged$seasonality)
seasonality_recode <- recode(seasonality_factor, "1" = "Rainy", "2" = "Dry")
seasonality_char <- as.character(seasonality_recode)
combined_df_merged$seasonality <- seasonality_char

combined_df_merged$province <- ifelse(
  is.na(combined_df_merged$seasonality),
  combined_df_merged$province,
  paste0(combined_df_merged$province, "_", combined_df_merged$seasonality)
)

# # JUST IN CASE... recalculate n.alleles for each locus of each sample
# combined_df_merged <- combined_df_merged %>%
#     group_by(NIDA2, locus) %>%
#     mutate(n.alleles = n_distinct(pseudo_cigar))

#recount n.alleles after filtering
combined_df_merged <- combined_df_merged %>%
  group_by(NIDA2, locus) %>%
  mutate(n.alleles = n_distinct(pseudo_cigar))

combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar) #rename alleles FOR DOWNSTREAM PROCESSING

combined_df_merged <- as.data.frame(combined_df_merged)

#remove Dry 
combined_df_merged <- combined_df_merged %>%
  filter(!grepl("Dry", province))


# FILTERING RESULTS
SS <- length(unique(combined_df_merged$NIDA2))
cat("Final sample size is", as.character(SS))

LC <- length(unique(combined_df_merged$locus))
cat("Final Diversity loci count is", as.character(LC))

sample_size_provinces <- combined_df_merged %>%
  group_by(province) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_provinces

sample_size_regions <- combined_df_merged %>%
  group_by(year, region) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_regions[,2:3]


#save allele_data_list
saveRDS(combined_df_merged, paste0("combined_df_merged_", SAMPLING, "_only.RDS"))
