
library(vegan)
library(dplyr)
library(RColorBrewer)

# assess suficiency of sampling

#input: combined_df_merged_SAMPLING_only.RDS

#output: accumulation_curves_provinces/regions_SAMPLING.pdf

######################################################

SAMPLING <- 2022 # 2021 or 2022

combined_df_merged <- readRDS(paste0("combined_df_merged_", SAMPLING, "_only.RDS")) 

combined_df_merged<- combined_df_merged[combined_df_merged$province %in% c("Maputo_Rainy", "Maputo_Dry", "Manica_Rainy", "Manica_Dry"),]

#######################################################
# 5.- SUFFICIENCY OF SAMPLINGS 

# ACCUMULATION CURVES (read this https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4885658/; ver genotype_curve de paquete poppr?)
#######################################################

# PROVINCE

# Initialize a list to store the rarefaction curves for each year
accum_curves_sub_SAMPLING <- list()

# Get unique years and provinces
unique_provinces <- unique(combined_df_merged$province)

# Iterate over each province
for (province in unique_provinces) {
  
  print(province)
  
  # Subsetting the data for SAMPLING
  sub_sub_SAMPLING <- combined_df_merged[combined_df_merged$year == SAMPLING & combined_df_merged$province == province, ]
  
  # Check if there are unique NIDA2s for SAMPLING
  if (length(unique(sub_sub_SAMPLING$NIDA2)) > 0) {
    # Initialize a list to store unique alleles for each NIDA2 for SAMPLING
    unique_alleles_sub_SAMPLING <- list()
    
    # Iterate over each unique NIDA2 for SAMPLING
    for (nida in unique(sub_sub_SAMPLING$NIDA2)) {
      subset_data <- sub_sub_SAMPLING[sub_sub_SAMPLING$NIDA2 == nida, ]
      unique_alleles_it <- unique(subset_data$allele)
      unique_alleles_sub_SAMPLING[[as.character(nida)]] <- unique_alleles_it
    }
    
    # Get unique alleles across all elements for SAMPLING
    all_unique_alleles_sub_SAMPLING <- unique(unlist(unique_alleles_sub_SAMPLING))
    
    # Create a matrix to store presence/absence of unique alleles for each element for SAMPLING
    presence_matrix_sub_SAMPLING <- sapply(unique_alleles_sub_SAMPLING, function(x) {
      as.integer(all_unique_alleles_sub_SAMPLING %in% x)
    })
    
    # Convert the matrix to a dataframe for SAMPLING
    presence_df_sub_SAMPLING <- as.data.frame(presence_matrix_sub_SAMPLING)
    presence_df_sub_SAMPLING <- t(presence_df_sub_SAMPLING)
    rownames(presence_df_sub_SAMPLING) <- names(unique_alleles_sub_SAMPLING)
    colnames(presence_df_sub_SAMPLING) <- all_unique_alleles_sub_SAMPLING
    
    # CALCULATE CURVE for SAMPLING
    accum_curve_sub_SAMPLING <- specaccum(presence_df_sub_SAMPLING, 'random', permutations = 100)
    accum_curves_sub_SAMPLING[[province]] <- accum_curve_sub_SAMPLING
  }
}

# Select 9 colors from the Paired palette
colors <- brewer.pal(11, "Paired")

pdf(paste0("accumulation_curves_provinces_", SAMPLING, ".pdf"), width = 12, height = 8)

max_samples <- max(sapply(accum_curves_sub_SAMPLING, function(x) length(x$sites)))
max_samples<- round(max_samples+(max_samples*0.05),0)

max_richness <- max(sapply(accum_curves_sub_SAMPLING, function(x) max(x$richness)))
max_richness<- round(max_richness+(max_richness*0.05),0)

# Plot the curves for SAMPLING
plot(accum_curves_sub_SAMPLING[[1]], col = colors[1], xlab = "Samples", ylab ="Alleles", main = "Allele Accumulation Curves per Province", xlim = c(0,max_samples), ylim = c(0,max_richness))
for (i in 2:length(accum_curves_sub_SAMPLING)) {
  lines(accum_curves_sub_SAMPLING[[i]], col = colors[i], lw = 1.5)
}
legend(x = max_samples/1.15, y = max_richness/3.7, legend = names(accum_curves_sub_SAMPLING), fill = colors, x.intersp = 0.7, y.intersp = 0.7)

dev.off()
#conclusion....