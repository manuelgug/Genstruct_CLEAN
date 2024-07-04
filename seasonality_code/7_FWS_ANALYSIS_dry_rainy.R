
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

SAMPLING <- 2022 # 2021 or 2022

combined_df_merged <- readRDS(paste0("combined_df_merged_", SAMPLING, "_only.RDS")) 
combined_df_merged<- combined_df_merged[combined_df_merged$province %in% c("Maputo_Rainy", "Maputo_Dry", "Manica_Rainy", "Manica_Dry"),]


#dry season is removed in the coi file, so no labels for dry
provinces <- c("Manica_Dry", "Manica_Rainy", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
province_colors <- c(Niassa = "firebrick4", Cabo_Delgado = "red", Nampula = "indianred1", Zambezia = "darkgreen", Tete = "forestgreen", Manica_Rainy = "springgreen2", Sofala  = "chartreuse", Inhambane = "cornflowerblue", Maputo_Rainy = "turquoise")



#######################################################
# 9.- He and Fws results 
#######################################################

# 1) calculate heterozygosity of the population (He); pop = province, region
#import everything into lists
rds_files <- list.files(path = ".", pattern = "\\MOIRE-RESULTS_FOR_ALLELE_FREQS.RDS$", full.names = TRUE)

# #remove dry
# if (SAMPLING == 2022){
#   
#   rds_files <- rds_files[!grepl("2021|Dry", rds_files)]
#   
# }else if (SAMPLING == 2021){
#   
#   rds_files <- rds_files[!grepl("2022", rds_files)]
#   
# }

rds_files <- rds_files[!rds_files %in% paste0("all_samples_complete_filtered_MOIRE-RESULTS_", SAMPLING, "_only_FOR_MOI.RDS")] 

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
processed_He_results$population <- gsub(paste0("_", SAMPLING, "_MOIRE-RESULTS_FOR_ALLELE_FREQS"), "", processed_He_results$population)

if (SAMPLING == 2021){
  
  processed_He_results <- processed_He_results[!(processed_He_results$population %in% c("Maputo", "Manica")), ] # remove MAPUTO because N too small
  
}

######
# keep amplicons with high He:

he_amps <- processed_He_results %>%
  group_by(locus) %>%
  summarize(mean = mean(post_stat_mean)) %>%
  arrange(desc(mean))

#keep top n% amplicons with highest He (or use a numeric threshold like > 0.1 of He or whatever... test.)
q95_amps<- round(length(unique(he_amps$locus))*0.95)
he_amps_q95 <- he_amps[1:q95_amps,]

# FILTER
processed_He_results <- processed_He_results[processed_He_results$locus %in% he_amps_q95$locus,]

library(stringr)
processed_He_results <- processed_He_results %>%
  mutate(geo = ifelse(str_detect(population, "North|South|Centre"), "region", "province"))




# 2) calculate heterozygosity of the individual (Hw): 𝐻W = 1 − (n𝑖(1/n𝑖)**2) 
heterozygosity_data <- combined_df_merged %>%
  group_by(NIDA2, locus) %>%
  mutate(Hw = 1 - (n.alleles * (1/n.alleles)^2))


# 3) calculate 1-Fws: 1-Fws = Hw/He
#merge He from provinces
merged_data <- heterozygosity_data %>%
  left_join(processed_He_results, by = c("locus" = "locus")) %>%
  filter(population == province) %>%
  mutate(He_province = ifelse(is.na(post_stat_mean), NA, post_stat_mean)) %>%
  select(NIDA2, locus, He_province)

heterozygosity_data <- heterozygosity_data %>%
  left_join(merged_data, by = c("NIDA2", "locus"))

heterozygosity_data <- distinct(heterozygosity_data)

#merge He from regions
merged_data <- heterozygosity_data %>%
  left_join(processed_He_results, by = c("locus" = "locus")) %>%
  filter(population == region) %>%
  mutate(He_region = ifelse(is.na(post_stat_mean), NA, post_stat_mean)) %>%
  select(NIDA2, locus, He_region)

heterozygosity_data <- heterozygosity_data %>%
  left_join(merged_data, by = c("NIDA2", "locus"))

heterozygosity_data <- distinct(heterozygosity_data)

#calculate 1-Fws for province and region as populations, for each year
heterozygosity_data$fws_province <- heterozygosity_data$Hw/heterozygosity_data$He_province
heterozygosity_data$fws_region <- heterozygosity_data$Hw/heterozygosity_data$He_region

#write.csv(heterozygosity_data, "He_and_Fws_2022.csv", row.names = F)

# Columns to keep
columns_to_keep <- c("NIDA2", "locus", "province", "region", "n.alleles",
                     "Hw", "He_province", "He_region", "fws_province", "fws_region")
# Filter columns
heterozygosity_data_filtered <- heterozygosity_data %>%
  select(all_of(columns_to_keep))

# Keep unique rows
heterozygosity_data_filtered <- distinct(heterozygosity_data_filtered)

# convert NA to 0 (monoallelic loci that doesn't have heterozygosity)
heterozygosity_data_filtered[is.na(heterozygosity_data_filtered)] <- 0

#sanity check
if ((length(unique(heterozygosity_data_filtered$NIDA2)) == length(unique(combined_df_merged$NIDA2))) & 
    (length(unique(heterozygosity_data_filtered$locus)) == length(unique(combined_df_merged$locus))) & 
    (length(unique(heterozygosity_data_filtered$province)) == length(unique(combined_df_merged$province))) & 
    (length(unique(heterozygosity_data_filtered$region)) == length(unique(combined_df_merged$region)))) {
  print("All looks good.")
}else{
  print("grab a coffee.")
}

#WHY IS IT NOT BETWEEN 0 AND 1??? (should it be?)
mean_Fws_per_individual<- heterozygosity_data_filtered %>%
  group_by(NIDA2, province) %>%
  summarize(mean_indiv_fws_province = mean(fws_province))

write.csv(mean_Fws_per_individual, paste0("mean_Fws_per_individual_dry_rainy_", SAMPLING, ".csv"), row.names = F)



mean_Fws_per_individual$province <- factor(mean_Fws_per_individual$province, levels = provinces)
mean_Fws_per_individual <- mean_Fws_per_individual%>%
  separate(province, into = c("province", "season"), sep = "_")

#kruskal
mean_Fws_per_individual %>%
  group_by(province) %>%
  summarise(kruskal_p_value = kruskal.test(mean_indiv_fws_province ~ season)$p.value)


mean_Fws_per_individual$season <- factor(mean_Fws_per_individual$season, levels = c("Rainy", "Dry"))
mean_Fws_per_individual$province <- factor(mean_Fws_per_individual$province, levels = c("Maputo", "Manica"))

# Retrieve the default ggplot2 color palette
default_colors <- scales::hue_pal()(2)

# Invert the order of the colors
inverted_colors <- default_colors

# Define the colors for the seasons, matching the order of levels
season_colors <- setNames(inverted_colors, c("Dry", "Rainy"))

prov_fws <- ggplot(mean_Fws_per_individual, aes(x = province, y = mean_indiv_fws_province, fill = season)) +
  geom_violin(width = 0.5, aes(color = season), alpha = 0.4) +
  geom_boxplot(width = 0.5, aes(color = season), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  labs(title = "",
       x = "",
       y = "Genome-wide 1-Fws")+
  scale_fill_manual(values = season_colors)+
  scale_color_manual(values = season_colors)
  
ggsave(paste0("province_fws_", SAMPLING,".png"), prov_fws, width = 8, height = 6, bg = "white")

