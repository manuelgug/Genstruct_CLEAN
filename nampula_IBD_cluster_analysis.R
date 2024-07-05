library(haven)
library(dcifer)
library(dplyr)
library(ggplot2)
library(reshape2)

merged_df <- readRDS("decifer_results_TABLE_2022.RDS")
combined_df_merged <- readRDS("combined_df_merged_2022_only.RDS")

#add VOC
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhps_doub_95_b")], by = "NIDA2")
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhfr_tr_95_b")], by = "NIDA2")
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhps_436_b")], by = "NIDA2")
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhps_581_b")], by = "NIDA2")


combined_df_merged$VOC <- ifelse(combined_df_merged$dhps_doub_95_b == 1 & combined_df_merged$dhfr_tr_95_b == 1, "dhps_d_dhfr_tr",
                                 ifelse(combined_df_merged$dhfr_tr_95_b == 1 & combined_df_merged$dhps_doub_95_b == 0, "dhfr_tr",
                                        ifelse(combined_df_merged$dhfr_tr_95_b == 0 & combined_df_merged$dhps_doub_95_b == 1, "dhps_d", "WT")))

#no_genotype
combined_df_merged$VOC <- ifelse(is.na(combined_df_merged$VOC), "no_genotype", combined_df_merged$VOC)

#436 581 genotypes
# Replace 0 with WT, 1 with mut, 2 with mix, and NA with no_Genotype in dhps_436
combined_df_merged$dhps_436 <- ifelse(combined_df_merged$dhps_436 == 0, "WT",
                                      ifelse(combined_df_merged$dhps_436 == 1, "mut",
                                             ifelse(combined_df_merged$dhps_436 == 2, "mix", "no_genotype")))

# Replace 0 with WT, 1 with mut, 2 with mix, and NA with no_Genotype in dhps_581
combined_df_merged$dhps_581 <- ifelse(combined_df_merged$dhps_581 == 0, "WT",
                                      ifelse(combined_df_merged$dhps_581 == 1, "mut",
                                             ifelse(combined_df_merged$dhps_581 == 2, "mix", "no_genotype")))

combined_df_merged$VOC_436_581 <- paste0("dhps_436", "-", combined_df_merged$dhps_436, "_", "581", "-", combined_df_merged$dhps_581)

combined_df_merged$VOC_436_581 <- gsub("dhps_436-NA_581-NA", "no_genotype", combined_df_merged$VOC_436_581 )
combined_df_merged$VOC_436_581 <- gsub("dhps_436-WT_581-WT", "WT", combined_df_merged$VOC_436_581 )

unique(combined_df_merged$VOC_436_581)

#remove NA rows
merged_df <- merged_df[complete.cases(merged_df),]

#keep significant cases
merged_df_signif <- merged_df[merged_df$p_value < 0.01,]
merged_df_signif <- merged_df[merged_df$estimate > 0.25,]

#merge with geo locations for each pair of samples compared
merged_df_signif_geo <- merge(merged_df_signif, combined_df_merged[, c("NIDA2", "province", "region", "VOC_436_581")], by.x = "sample1", by.y = "NIDA2")
colnames(merged_df_signif_geo)[7:9] <- c("province_s1", "region_s1", "VOC_s1")
merged_df_signif_geo<- distinct(merged_df_signif_geo)

merged_df_signif_geo <- merge(merged_df_signif_geo, combined_df_merged[, c("NIDA2", "province", "region", "VOC_436_581")], by.x = "sample2", by.y = "NIDA2")
colnames(merged_df_signif_geo)[10:12] <- c("province_s2", "region_s2", "VOC_s2")
merged_df_signif_geo<- distinct(merged_df_signif_geo)




# now analyze the clusters 
colnames(merged_df_signif_geo)

NampvsNamp_IBD <- merged_df_signif_geo[merged_df_signif_geo$province_s1 == "Nampula" & merged_df_signif_geo$province_s2 == "Nampula",]
NampvsNamp_IBD$sample2 <- as.character(NampvsNamp_IBD$sample2)
NampvsNamp_IBD$sample1 <- as.character(NampvsNamp_IBD$sample1)
NampvsNamp_IBD

db <- read_dta('DrugRes_2021_2022_DB_ALLDATA_1Jan2024.dta', )
db_subset <- db[c("NIDA2", "district")]

# Merge the first data frame based on sample1 and NIDA2
merged_df <- merge(NampvsNamp_IBD, db_subset, by.x = "sample1", by.y = "NIDA2", all.x = TRUE)

# Rename the district column to district1
names(merged_df)[names(merged_df) == "district"] <- "district1"

# Merge the resulting data frame based on sample2 and NIDA2
final_merged_df <- merge(merged_df, db_subset, by.x = "sample2", by.y = "NIDA2", all.x = TRUE)

# Rename the district column to district2
names(final_merged_df)[names(final_merged_df) == "district"] <- "district2"


NampvsNamp_IBD <- final_merged_df

NampvsNamp_IBD <- NampvsNamp_IBD[order(-NampvsNamp_IBD$estimate), ]

NampvsNamp_IBD <- NampvsNamp_IBD %>%
  select(sample1, sample2, everything())



#add coi data
coi <- read.csv("coi_for_db_2022.csv")
coi$NIDA2 <- as.character(coi$NIDA2)

coi_subset <- coi[c("NIDA2", "post_effective_coi_med")]

# Merge the first data frame based on sample1 and NIDA2
merged_df <- merge(NampvsNamp_IBD, coi_subset, by.x = "sample1", by.y = "NIDA2", all.x = TRUE)

# Rename the district column to district1
names(merged_df)[names(merged_df) == "post_effective_coi_med"] <- "ecoi1"

# Merge the resulting data frame based on sample2 and NIDA2
final_merged_df <- merge(merged_df, coi_subset, by.x = "sample2", by.y = "NIDA2", all.x = TRUE)

# Rename the district column to district2
names(final_merged_df)[names(final_merged_df) == "post_effective_coi_med"] <- "ecoi2"


NampvsNamp_IBD <- final_merged_df

NampvsNamp_IBD <- NampvsNamp_IBD[order(-NampvsNamp_IBD$estimate), ]

NampvsNamp_IBD <- NampvsNamp_IBD %>%
  select(sample1, sample2, everything())


write.csv(NampvsNamp_IBD, "NampvsNamp_IBD.csv",row.names = F)


clone_like <- NampvsNamp_IBD[NampvsNamp_IBD$estimate > 0.9,]

as.data.frame(table(paste0(clone_like$district1, "_", clone_like$district2)))


# Combine sample1 and sample2 into one column
combined_samples <- clone_like %>%
  mutate(combined = paste(sample1, sample2, sep = ", ")) %>%
  select(-sample1, -sample2)

library(tidyr)
# Count unique strings in sample1 and sample2 for each combination of district1 and district2
unique_counts <- clone_like %>%
  gather(key = "sample_type", value = "sample", sample1, sample2) %>%
  group_by(district1, district2) %>%
  summarize(unique_count = n_distinct(sample)) %>%
  ungroup()
