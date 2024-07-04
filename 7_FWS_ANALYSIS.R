
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

SAMPLING <- 2021 # 2021 or 2022

combined_df_merged <- readRDS(paste0("combined_df_merged_", SAMPLING, "_only.RDS")) 
combined_df_merged <- combined_df_merged[!(combined_df_merged$province %in% c("Maputo_Dry", "Manica_Dry")), ] # remove dry



#set labels
if (SAMPLING == 2022){
  
  #dry season is removed in the coi file, so no labels for dry
  provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Rainy") #ordered from north to south
  province_colors <- c(Niassa = "firebrick4", Cabo_Delgado = "red", Nampula = "indianred1", Zambezia = "darkgreen", Tete = "forestgreen", Manica_Rainy = "springgreen2", Sofala  = "chartreuse", Inhambane = "cornflowerblue", Maputo_Rainy = "turquoise")
  
  
} else if (SAMPLING == 2021){
  
  combined_df_merged <- combined_df_merged[!(combined_df_merged$province %in% c("Maputo", "Manica")), ] # remove dry
  
  provinces <- c("Niassa", "Nampula", "Zambezia", "Inhambane") #ordered from north to south "Maputo" and "Manica" out because small N
  province_colors <- c(Niassa = "firebrick4", Nampula = "indianred1", Zambezia = "darkgreen",  Inhambane = "cornflowerblue") #Maputo = "deepskyblue", Manica = "green",
  
}else{
  
  print(paste0("No samples for year ", SAMPLING))
  
}


#######################################################
# 9.- He and Fws results 
#######################################################

# 1) calculate heterozygosity of the population (He); pop = province, region
#import everything into lists
rds_files <- list.files(path = ".", pattern = "\\MOIRE-RESULTS_FOR_ALLELE_FREQS.RDS$", full.names = TRUE)

#remove dry
if (SAMPLING == 2022){
  
  rds_files <- rds_files[!grepl("2021|Dry", rds_files)]
  
}else if (SAMPLING == 2021){
  
  rds_files <- rds_files[!grepl("2022", rds_files)]
  
}

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
  group_by(NIDA2, province, region) %>%
  summarize(mean_indiv_fws_province = mean(fws_province),
            mean_indiv_fws_region = mean(fws_region))

write.csv(mean_Fws_per_individual, paste0("mean_Fws_per_individual_", SAMPLING, ".csv"), row.names = F)


regions <- c("North", "Centre", "South")

mean_Fws_per_individual$province <- factor(mean_Fws_per_individual$province, levels = provinces)
mean_Fws_per_individual$region <- factor(mean_Fws_per_individual$region, levels = regions)


pairwise_province_fws <- pairwise.wilcox.test(mean_Fws_per_individual$mean_indiv_fws_province, 
                                              mean_Fws_per_individual$province, p.adjust.method = "bonferroni")

pairwise_province_fws <- melt(pairwise_province_fws[[3]])
signif_p.pairwise_province_fws<- pairwise_province_fws[pairwise_province_fws$value <0.05 & !is.na(pairwise_province_fws$value),]

combos_pairwise_province_fws_province <- lapply(1:nrow(signif_p.pairwise_province_fws), function(i) {
  as.character(c(signif_p.pairwise_province_fws[i, "Var1"], 
                 signif_p.pairwise_province_fws[i, "Var2"]))
})

if (SAMPLING == 2021){
  
  mean_Fws_per_individual_nomaputo <- mean_Fws_per_individual[!mean_Fws_per_individual$province == "Maputo",] #MAPUTO'S FWS IN 2021 IS HUGE BEACUSE SAMPLE SIZE IS 4. EXCLUDE

  prov_fws <- ggplot(mean_Fws_per_individual_nomaputo, aes(x = province, y = mean_indiv_fws_province, fill = province)) +
    geom_violin(width = 1, aes(color = region), alpha = 0.4) +
    geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11), 
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_discrete(name = "Region") +  # Customize legend title
    #ggtitle("Province Connectivity") +
    labs(x = "Province", y = "Genome-wide 1-Fws") +
    guides(color = FALSE, fill = FALSE) +
    scale_fill_manual(values = province_colors) #+
    # stat_compare_means(comparisons = combos_pairwise_province_fws_province, aes(label = after_stat(p.signif)),
    #                    method = "wilcox.test")
  
  
  # Prepare significance labels
  signif_p.pairwise_province_fws <- signif_p.pairwise_province_fws %>%
    mutate(label = case_when(
      value < 0.001 ~ "***",
      value < 0.01 ~ "**",
      value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  signif_p.pairwise_province_fws <- signif_p.pairwise_province_fws[signif_p.pairwise_province_fws$value < 0.05,]
  
  
  # Add significance annotations only for significant results
  for (i in 1:nrow(signif_p.pairwise_province_fws)) {
    if (signif_p.pairwise_province_fws$label[i] != "") {
      prov_fws <- prov_fws + geom_signif(
        comparisons = list(c(as.character(signif_p.pairwise_province_fws$Var1)[i], as.character(signif_p.pairwise_province_fws$Var2)[i])),
        annotations = signif_p.pairwise_province_fws$label[i],
        map_signif_level = TRUE,
        y_position = max(mean_Fws_per_individual_nodry$mean_indiv_fws_region) + 0.1 * (i + 0.1)
      )
    }
  }
  
  prov_fws
  
  
  ggsave(paste0("province_fws_", SAMPLING,".png"), prov_fws, width = 8, height = 6, bg = "white")
  
} else{
  
  prov_fws <- ggplot(mean_Fws_per_individual, aes(x = province, y = mean_indiv_fws_province, fill = province)) +
    geom_violin(width = 1, aes(color = region), alpha = 0.4) +
    geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11), 
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_discrete(name = "Region") +  # Customize legend title
    #ggtitle("Province Connectivity") +
    labs(x = "Province", y = "Genome-wide 1-Fws") +
    guides(color = FALSE, fill = FALSE) +
    scale_fill_manual(values = province_colors) #+
    # stat_compare_means(comparisons = combos_pairwise_province_fws_province, aes(label = after_stat(p.signif)),
    #                    method = "wilcox.test")
  
  
  # Prepare significance labels
  signif_p.pairwise_province_fws <- signif_p.pairwise_province_fws %>%
    mutate(label = case_when(
      value < 0.001 ~ "***",
      value < 0.01 ~ "**",
      value < 0.05 ~ "*",
      TRUE ~ ""
    ))
  signif_p.pairwise_province_fws <- signif_p.pairwise_province_fws[signif_p.pairwise_province_fws$value < 0.05,]
  
  
  # Add significance annotations only for significant results
  for (i in 1:nrow(signif_p.pairwise_province_fws)) {
    if (signif_p.pairwise_province_fws$label[i] != "") {
      prov_fws <- prov_fws + geom_signif(
        comparisons = list(c(as.character(signif_p.pairwise_province_fws$Var1)[i], as.character(signif_p.pairwise_province_fws$Var2)[i])),
        annotations = signif_p.pairwise_province_fws$label[i],
        map_signif_level = TRUE,
        y_position = max(mean_Fws_per_individual_nodry$mean_indiv_fws_region) + 0.1 * (i + 0.1)
      )
    }
  }
  
  prov_fws
  
  ggsave(paste0("province_fws_", SAMPLING,".png"), prov_fws, width = 8, height = 6, bg = "white")
  
}




mean_Fws_per_individual_nodry <- mean_Fws_per_individual %>%
  filter(!grepl("Dry", province))

pairwise_region_fws <- pairwise.wilcox.test(mean_Fws_per_individual_nodry$mean_indiv_fws_region, 
                                            mean_Fws_per_individual_nodry$region, p.adjust.method = "bonferroni")

pairwise_region_fws <- melt(pairwise_region_fws[[3]])
signif_p.pairwise_region_fws<- pairwise_region_fws[pairwise_region_fws$value <0.05 & !is.na(pairwise_region_fws$value),]

combos_pairwise_region_fws_province <- lapply(1:nrow(signif_p.pairwise_region_fws), function(i) {
  as.character(c(signif_p.pairwise_region_fws[i, "Var1"], 
                 signif_p.pairwise_region_fws[i, "Var2"]))
})

reg_fws <- ggplot(mean_Fws_per_individual_nodry, aes(x = region, y = mean_indiv_fws_region, fill = region)) +
  geom_violin(width = 1, aes(color = region), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Genome-wide 1-Fws") +
  guides(fill = FALSE, color = FALSE) #+
  #stat_compare_means(comparisons = combos_pairwise_region_fws_province, aes(label = after_stat(p.signif)),
  #                  method = "wilcox.test")


# Prepare significance labels
signif_p.pairwise_region_fws <- signif_p.pairwise_region_fws %>%
  mutate(label = case_when(
    value < 0.001 ~ "***",
    value < 0.01 ~ "**",
    value < 0.05 ~ "*",
    TRUE ~ ""
  ))
signif_p.pairwise_region_fws <- signif_p.pairwise_region_fws[signif_p.pairwise_region_fws$value < 0.05,]


# Add significance annotations only for significant results
for (i in 1:nrow(signif_p.pairwise_region_fws)) {
  if (signif_p.pairwise_region_fws$label[i] != "") {
    reg_fws <- reg_fws + geom_signif(
      comparisons = list(c(as.character(signif_p.pairwise_region_fws$Var1)[i], as.character(signif_p.pairwise_region_fws$Var2)[i])),
      annotations = signif_p.pairwise_region_fws$label[i],
      map_signif_level = TRUE,
      y_position = max(mean_Fws_per_individual_nodry$mean_indiv_fws_region) + 0.1 * (i + 0.1)
    )
  }
}

ggsave(paste0("region_fws_", SAMPLING,".png"), reg_fws, width = 8, height = 6, bg = "white")
