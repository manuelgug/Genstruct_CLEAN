
library(dplyr)
library(cowplot)
library(ggplot2)
library(vegan)
library(ape)    
library(geosphere)
library(Rtsne)
library(reshape2)
library(haven)

SAMPLING <- 2022 # 2021 or 2022

combined_df_merged <- readRDS(paste0("combined_df_merged_", SAMPLING, "_only.RDS")) 
combined_df_merged <- combined_df_merged[!(combined_df_merged$province %in% c("Maputo_Dry", "Manica_Dry")), ] # remove dry

db <- read_dta('DrugRes_2021_2022_DB_ALLDATA_1Jan2024.dta')


#set labels
if (SAMPLING == 2022){
  
  #dry season is removed in the coi file, so no labels for dry
  provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Rainy") #ordered from north to south
  province_colors <- c(Niassa = "firebrick4", Cabo_Delgado = "red", Nampula = "indianred1", Zambezia = "darkgreen", Tete = "forestgreen", Manica_Rainy = "springgreen2", Sofala  = "chartreuse", Inhambane = "cornflowerblue", Maputo_Rainy = "turquoise")
  
  
} else if (SAMPLING == 2021){
  
  provinces <- c("Niassa", "Nampula", "Zambezia", "Manica", "Inhambane", "Maputo") #ordered from north to south
  province_colors <- c(Niassa = "firebrick4", Nampula = "indianred1", Zambezia = "darkgreen", Manica = "green", Inhambane = "cornflowerblue", Maputo = "deepskyblue")
  
}else{
  
  print(paste0("No samples for year ", SAMPLING))
  
}





########################
# 11. Ordination, population structure
########################

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


#input for multivariate analyses
raref_input <- as.data.frame(cbind(NIDA2 = combined_df_merged$NIDA2, 
                                   year = combined_df_merged$year, 
                                   province = combined_df_merged$province,
                                   region = combined_df_merged$region,
                                   locus = combined_df_merged$locus,
                                   n.alleles = combined_df_merged$n.alleles,
                                   norm.reads.locus = combined_df_merged$norm.reads.locus,
                                   allele = paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar),
                                   run_id = combined_df_merged$run_id,
                                   VOC = combined_df_merged$VOC,
                                   VOC_436_581 = combined_df_merged$VOC_436_581))


#remove Dry
raref_input <- raref_input %>%
  filter(!grepl("Dry", province))


# Melt the data frame to convert it from wide to long format
melted <- melt(raref_input, id.vars = c("NIDA2", "norm.reads.locus"), measure.vars = "allele")
melted<-melted[,-3]

library(tidyr)

rearranged <- melted %>%
  pivot_wider(
    names_from = value,
    values_from = norm.reads.locus
  )

# format
rearranged <- as.data.frame(rearranged)
rownames(rearranged) <- rearranged$NIDA2
rearranged <- rearranged[, -1]
rearranged <- replace(rearranged, is.na(rearranged), 0)

#pca labels:
NIDA2 <-data.frame(NIDA2 = rownames(rearranged))
pca_labels<- combined_df_merged %>% distinct(NIDA2, year, province, region, VOC_436_581)

pca_labels <- pca_labels %>%
  filter(!grepl("Dry", province))

if (all(NIDA2$NIDA2 == pca_labels$NIDA2)){
  print("Order of categorical variables is ok.")
}else{
  "grab a coffee."
}

# format freq df
rearranged <- rearranged %>%
  mutate_all(as.numeric)

zero_cols <- sapply(rearranged, function(x) all(x == 0))
rearranged_filtered <- rearranged[, !zero_cols]

# Replace values greater than 0 with 1: MAKE IT PRESENCE/ABSENCE
rearranged_pres_abs <- rearranged %>%
  mutate_all(~ ifelse(. > 0, 1, .))

rearranged_pres_abs <- rearranged_pres_abs %>%
  mutate_all(as.numeric)



#PCoA

if (SAMPLING == 2022){
  
  shapes <- c(13, 12, 11, 1, 19)
  shapes2 <- c(13, 1, 19)
  
}else if (SAMPLING == 2021){
  
  shapes <-  c(13, 12, 11, 10, 9, 19)
  shapes2 <- c(13, 10, 11, 19)
}

# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(rearranged, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract the principal coordinate scores
pcs <- as.data.frame(pcoa_result$vectors)

# Combine principal coordinate scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 = pca_labels$VOC_436_581)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)

regions <- c("North", "Centre", "South")
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# # Plot PCoA
variance_explained <- round(pcoa_result$values$Eigenvalues / sum(pcoa_result$values$Eigenvalues) * 100, 2)
variance_explained_axis1 <- variance_explained[1]
variance_explained_axis2 <- variance_explained[2]


# Plot PCoA with variance explained in title
af_pcoa <- ggplot(pcs_with_labels, aes(x = Axis.1, y = Axis.2, color = province, shape = VOC_436_581)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()+
  theme(
    legend.position = "bottom", 
    legend.box = "vertical",  
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.byrow = TRUE,  # Place legends by row
    axis.title.x = element_text(size = 14),  # Increase x-axis label size
    axis.title.y = element_text(size = 14)) +  # Increase y-axis label size
  guides(fill = FALSE, color = FALSE, shape = FALSE)+
  scale_color_manual(values = province_colors)+
  scale_shape_manual(values = shapes)

af_pcoa

#PCoA presence/absence

# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(rearranged_pres_abs, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract the principal coordinate scores
pcs <- as.data.frame(pcoa_result$vectors)

# Combine principal coordinate scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 = pca_labels$VOC_436_581)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)


# # Plot PCoA
variance_explained <- round(pcoa_result$values$Eigenvalues / sum(pcoa_result$values$Eigenvalues) * 100, 2)
variance_explained_axis1 <- variance_explained[1]
variance_explained_axis2 <- variance_explained[2]

# Plot PCoA with variance explained in title
pa_pcoa <- ggplot(pcs_with_labels, aes(x = Axis.1, y = Axis.2, color = province, shape = VOC_436_581)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()+
  theme(
    legend.position = "bottom", 
    legend.box = "vertical",  
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.byrow = F,  # Place legends by row
    axis.title.x = element_text(size = 14),  # Increase x-axis label size
    axis.title.y = element_text(size = 14)   # Increase y-axis label size
  ) +
  scale_color_manual(values = province_colors)+
  scale_shape_manual(values = shapes)+
  guides(
    color = guide_legend(ncol = 3),  # Arrange color legend in 3 columns
    shape = guide_legend(ncol = 3)   # Arrange shape legend in 3 columns
  )

pa_pcoa


combined_plot_pcoa <- plot_grid(af_pcoa, pa_pcoa, ncol = 2)

ggsave(paste0("PCoA_regions_", SAMPLING, ".png"), combined_plot_pcoa, width = 16, height = 10, bg = "white")



# pcoa separated by mutations

#format 436 mutation labels
pcs_with_labels$VOC_436 <- ifelse(pcs_with_labels$VOC_436_581 == "dhps_436-mix_581-WT" | pcs_with_labels$VOC_436_581 == "dhps_436-mut_581-WT", 
                                  "dhps_436_mut/mix", pcs_with_labels$VOC_436_581)

pcs_with_labels$VOC_436 <- ifelse(pcs_with_labels$VOC_436 == "dhps_436-WT_581-mix", "WT", pcs_with_labels$VOC_436)

#format 581 mutations labels
pcs_with_labels$VOC_581 <- ifelse(pcs_with_labels$VOC_436_581 == "dhps_436-mix_581-WT" | pcs_with_labels$VOC_436_581 == "dhps_436-mut_581-WT", 
                                  "WT", pcs_with_labels$VOC_436_581)

pcs_with_labels$VOC_581 <- ifelse(pcs_with_labels$VOC_581 == "dhps_436-WT_581-mix", "dhps_581_mix", pcs_with_labels$VOC_581)


#PCoA 436
# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(rearranged, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract the principal coordinate scores
pcs <- as.data.frame(pcoa_result$vectors)

# Combine principal coordinate scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 = pca_labels$VOC_436_581)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# # Plot PCoA
variance_explained <- round(pcoa_result$values$Eigenvalues / sum(pcoa_result$values$Eigenvalues) * 100, 2)
variance_explained_axis1 <- variance_explained[1]
variance_explained_axis2 <- variance_explained[2]

#format 436 mutation labels
pcs_with_labels$VOC_436 <- ifelse(pcs_with_labels$VOC_436_581 == "dhps_436-mix_581-WT" | pcs_with_labels$VOC_436_581 == "dhps_436-mut_581-WT", 
                                  "dhps_436_mut/mix", pcs_with_labels$VOC_436_581)

pcs_with_labels$VOC_436 <- ifelse(pcs_with_labels$VOC_436 == "dhps_436-WT_581-mix", "WT", pcs_with_labels$VOC_436)

# Plot PCoA with variance explained in title
af_pcoa <- ggplot(pcs_with_labels, aes(x = Axis.1, y = Axis.2, color = province, shape = VOC_436)) +
  geom_point(size = 4, alpha = 0.5) +
  labs(title = "",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()+
  guides(fill = FALSE, color = FALSE, shape = FALSE)+
  scale_color_manual(values = province_colors)+
  scale_shape_manual(values = shapes2)

af_pcoa

#PCoA presence/absence

# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(rearranged_pres_abs, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract the principal coordinate scores
pcs <- as.data.frame(pcoa_result$vectors)

# Combine principal coordinate scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 = pca_labels$VOC_436_581)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# # Plot PCoA
variance_explained <- round(pcoa_result$values$Eigenvalues / sum(pcoa_result$values$Eigenvalues) * 100, 2)
variance_explained_axis1 <- variance_explained[1]
variance_explained_axis2 <- variance_explained[2]

#format 436 mutation labels
pcs_with_labels$VOC_436 <- ifelse(pcs_with_labels$VOC_436_581 == "dhps_436-mix_581-WT" | pcs_with_labels$VOC_436_581 == "dhps_436-mut_581-WT", 
                                  "dhps_436_mut/mix", pcs_with_labels$VOC_436_581)

pcs_with_labels$VOC_436 <- ifelse(pcs_with_labels$VOC_436 == "dhps_436-WT_581-mix", "WT", pcs_with_labels$VOC_436)

# Plot PCoA with variance explained in title
pa_pcoa <- ggplot(pcs_with_labels, aes(x = Axis.1, y = Axis.2, color = province, shape = VOC_436)) +
  geom_point(size = 4, alpha = 0.5) +
  labs(title = "",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()+
  scale_color_manual(values = province_colors)+
  scale_shape_manual(values = shapes2)

pa_pcoa

combined_plot_pcoa <- plot_grid(af_pcoa, pa_pcoa, ncol = 2)

ggsave(paste0("PCoA_regions_436_", SAMPLING, ".png"), combined_plot_pcoa, width = 16, height = 10, bg = "white")



#PCoA 581
# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(rearranged, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract the principal coordinate scores
pcs <- as.data.frame(pcoa_result$vectors)

# Combine principal coordinate scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 = pca_labels$VOC_436_581)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# # Plot PCoA
variance_explained <- round(pcoa_result$values$Eigenvalues / sum(pcoa_result$values$Eigenvalues) * 100, 2)
variance_explained_axis1 <- variance_explained[1]
variance_explained_axis2 <- variance_explained[2]

#format 581 mutations labels
pcs_with_labels$VOC_581 <- ifelse(pcs_with_labels$VOC_436_581 == "dhps_436-mix_581-WT" | pcs_with_labels$VOC_436_581 == "dhps_436-mut_581-WT", 
                                  "WT", pcs_with_labels$VOC_436_581)

pcs_with_labels$VOC_581 <- ifelse(pcs_with_labels$VOC_581 == "dhps_436-WT_581-mix", "dhps_581_mix", pcs_with_labels$VOC_581)

# Plot PCoA with variance explained in title
af_pcoa <- ggplot(pcs_with_labels, aes(x = Axis.1, y = Axis.2, color = province, shape = VOC_581)) +
  geom_point(size = 4, alpha = 0.5) +
  labs(title = "",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()+
  guides(fill = FALSE, color = FALSE, shape = FALSE)+
  scale_color_manual(values = province_colors)+
  scale_shape_manual(values = shapes2)

af_pcoa

#PCoA presence/absence

# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(rearranged_pres_abs, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract the principal coordinate scores
pcs <- as.data.frame(pcoa_result$vectors)

# Combine principal coordinate scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 = pca_labels$VOC_436_581)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# # Plot PCoA
variance_explained <- round(pcoa_result$values$Eigenvalues / sum(pcoa_result$values$Eigenvalues) * 100, 2)
variance_explained_axis1 <- variance_explained[1]
variance_explained_axis2 <- variance_explained[2]

#format 581 mutations labels
pcs_with_labels$VOC_581 <- ifelse(pcs_with_labels$VOC_436_581 == "dhps_436-mix_581-WT" | pcs_with_labels$VOC_436_581 == "dhps_436-mut_581-WT", 
                                  "WT", pcs_with_labels$VOC_436_581)

pcs_with_labels$VOC_581 <- ifelse(pcs_with_labels$VOC_581 == "dhps_436-WT_581-mix", "dhps_581_mix", pcs_with_labels$VOC_581)

# Plot PCoA with variance explained in title
pa_pcoa <- ggplot(pcs_with_labels, aes(x = Axis.1, y = Axis.2, color = province, shape = VOC_581)) +
  geom_point(size = 4, alpha = 0.5) +
  labs(title = "",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()+
  scale_color_manual(values = province_colors)+
  scale_shape_manual(values = shapes2)

pa_pcoa

combined_plot_pcoa <- plot_grid(af_pcoa, pa_pcoa, ncol = 2)

ggsave(paste0("PCoA_regions_581_", SAMPLING, ".png"), combined_plot_pcoa, width = 16, height = 10, bg = "white")



# POPULATION ALLELE FREQ TSNE

# 1) extract allele freqs
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
allele_freqs_list <- list()

for (file in rds_files) {
  print(file)
  file_name <- tools::file_path_sans_ext(basename(file))
  allele_freqs_list[[file_name]] <- readRDS(file)
}

# Loop through each element in allele_freqs_list
processed_allele_freq_results <- data.frame()

for (i in seq_along(allele_freqs_list)) {
  
  # Summarize He
  allele_freq_results <- moire::summarize_allele_freqs(allele_freqs_list[[i]])
  allele_freq_results$population <- names(allele_freqs_list[i])
  
  # Add the processed He results to the list
  processed_allele_freq_results <- rbind(processed_allele_freq_results, allele_freq_results)
}

#formatting categories
processed_allele_freq_results$population <- gsub(paste0("_", SAMPLING,"_MOIRE-RESULTS_FOR_ALLELE_FREQS"), "", processed_allele_freq_results$population)

library(stringr)
processed_allele_freq_results <- processed_allele_freq_results %>%
  mutate(geo = ifelse(str_detect(population, "North|South|Centre"), "region", "province"))


# Melt the data frame to convert it from wide to long format
melted <- melt(processed_allele_freq_results, id.vars = c("population", "post_allele_freqs_mean"), measure.vars = "allele")
melted<-melted[,-3]

library(tidyr)

rearranged_processed_allele_freq_results <- melted %>%
  pivot_wider(
    names_from = value,
    values_from = post_allele_freqs_mean
  )

# format
rearranged_processed_allele_freq_results <- as.data.frame(rearranged_processed_allele_freq_results)
rownames(rearranged_processed_allele_freq_results) <- rearranged_processed_allele_freq_results$population
rearranged_processed_allele_freq_results <- rearranged_processed_allele_freq_results[, -1]
rearranged_processed_allele_freq_results <- replace(rearranged_processed_allele_freq_results, is.na(rearranged_processed_allele_freq_results), 0)


# Find rows with "Centre", "North", or "South" in their names
region_columns <- grepl("Centre|North|South", rownames(rearranged_processed_allele_freq_results))
rearranged_processed_allele_freq_results_region <- rearranged_processed_allele_freq_results[region_columns, ]
rearranged_processed_allele_freq_results_province <- rearranged_processed_allele_freq_results[!region_columns, ]

library(stringr)
# Split row names by the last "_"
province_region_eq <- data.frame(
  metadata_province = c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Rainy"),
  metadata_region = c("North", "North", "North", "Centre", "Centre", "Centre", "Centre", "South", "South")
)

province_region_eq <- rbind(province_region_eq, c("Maputo", "South"), c("Manica", "Centre"))
metadata_province <- rownames(rearranged_processed_allele_freq_results_province)
metadata_province <- factor(metadata_province, levels = provinces)
metadata_mapping <- merge(data.frame(metadata_province), province_region_eq, by = "metadata_province", all.x = TRUE)
metadata_region <- metadata_mapping$metadata_region



# tsne of provinces
perplexity <- floor((length(metadata_province) - 1) / 3) #highest possible
set.seed(69)
tsne_result_freqs <- Rtsne(as.matrix(rearranged_processed_allele_freq_results_province), dims = 2, verbose = TRUE, check_duplicates = FALSE, pca_center = T, max_iter = 2e4, num_threads = 0, perplexity = perplexity)

# Convert t-SNE results to data frame
tsne_data_freqs <- as.data.frame(tsne_result_freqs$Y)

tsne_data_freqs<- cbind(tsne_data_freqs, metadata_region)

# Plot t-SNE of freqs
pop_allele_freq_tsne <- ggplot(tsne_data_freqs, aes(V1, V2, color = metadata_province, shape = metadata_region)) + # shape = factor(metadata_province$site)
  geom_point(size = 8, alpha = 0.7) +
  labs(title = "",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()+
  scale_color_manual(values = province_colors)

pop_allele_freq_tsne

ggsave(paste0("pop_allele_freqs_tsne", SAMPLING, ".png"), pop_allele_freq_tsne, width = 10, height = 8, bg = "white")






### MANTEL TEST

# Remove rows based on the indices and reorder rows
rows_to_remove <- grepl("Dry", rownames(rearranged_processed_allele_freq_results_province))
rearranged_processed_allele_freq_results_province_nodry <- rearranged_processed_allele_freq_results_province[!rows_to_remove, ]

#reorder rows
indices <- match(provinces, rownames(rearranged_processed_allele_freq_results_province_nodry))
rearranged_processed_allele_freq_results_province_nodry <- rearranged_processed_allele_freq_results_province_nodry[indices, ]

#calculate distances
bray_curtis_dist <- vegdist(rearranged_processed_allele_freq_results_province_nodry, method = "bray", diag = T, upper = T)
# manhattan_dist <- vegdist(rearranged_processed_allele_freq_results_province_nodry, method = "manhattan", diag = T, upper = T)
# euclidean_dist <- vegdist(rearranged_processed_allele_freq_results_province_nodry, method = "euclidean", diag = T, upper = T)
# gower_dist <- vegdist(rearranged_processed_allele_freq_results_province_nodry, method = "gower", diag = T, upper = T)

# not actual coordinates??
coordinates <- data.frame(
  longitude = c(36.4931, 39.6123, 38.3473, 36.8663, 33.5867, 33.898, 33.898, 34.8444, 35.3837, 32.5716, 32.5716),
  latitude = c(-12.8779, -11.6455, -15.1056, -16.2828, -16.1742, -19.0834, -19.0834, -19.1548, -23.865, -25.9692, -25.9692)
)

rownames(coordinates) <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Rainy", "Manica", "Sofala", "Inhambane", "Maputo_Rainy", "Maputo")

coordinates<- coordinates[rownames(coordinates) %in% provinces, ]

#calculate harvesine distance
geo_dist <- distm(coordinates, fun = distHaversine)
rownames(geo_dist) <- provinces
colnames(geo_dist) <- provinces


#mantel test fucntion to test multiple distance metrics
perform_mantel_test <- function(distance, geo_dist) {
  
  # Perform Mantel test
  mantel_result <- mantel(distance, geo_dist, method = "spearman", permutations = 10000)
  
  # Extract distance vectors
  distance_vec <- as.vector(as.matrix(distance))
  geo_vec <- as.vector(geo_dist)
  
  # Create data frame for plotting
  mat <- data.frame(deest = distance_vec, geo = geo_vec)
  
  # Filter out zero values for Bray-Curtis distance
  mat <- mat[mat$deest > 0, ]
  
  # Create scatter plot
  plot <- ggplot(mat, aes(y = deest, x = geo / 1000)) + 
    geom_point(size = 4, alpha = 0.75, colour = "black", shape = 21, aes(fill = geo / 1000)) + 
    geom_smooth(method = "lm", colour = "black", alpha = 0.2) + 
    labs(x = "Harvesine Distance", y = "Bray-Curtis Dissimilarity", fill = "Kilometers") + 
    theme(axis.text.x = element_text(colour = "black", size = 12), 
          axis.text.y = element_text(size = 11, colour = "black"), 
          axis.title = element_text(size = 14, colour = "black"), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "right",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11)) +
    scale_fill_continuous(high = "navy", low = "skyblue")
  
  # Return Mantel test result and the scatter plot
  return(list(mantel_result = mantel_result, plot = plot))
}

res <- perform_mantel_test(bray_curtis_dist, geo_dist)
# perform_mantel_test(manhattan_dist, geo_dist)
# perform_mantel_test(euclidean_dist, geo_dist)
# perform_mantel_test(gower_dist, geo_dist)

p <- res$plot

pval <- round(res$mantel_result$signif, 3)
mantelsR <- round(res$mantel_result$statistic,3)

p <- p + 
  annotate("text", x = Inf, y = -Inf, label = paste("Mantel's R:", mantelsR, "\n", "p-value:", pval), 
           hjust = 1.1, vjust = -0.5, size = 5, color = "black")

ggsave(paste0("mantel_pop_allele_Freqs", SAMPLING,".png"), p, width = 8, height = 6, bg = "white")
