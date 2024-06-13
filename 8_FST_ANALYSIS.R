
library(dplyr)
library(ggplot2)

SAMPLING <- 2022 # 2021 or 2022

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


if (SAMPLING == 2021){
  
  processed_He_results <- processed_He_results[!(processed_He_results$population %in% c("Maputo", "Manica")), ] # remove MAPUTO because N too small
  
}


#separate provinces and regions
processed_He_results_provinces <- processed_He_results[processed_He_results$geo == "province", ]
processed_He_results_provinces$pop <- processed_He_results_provinces$population
processed_He_results_regions<- processed_He_results[processed_He_results$geo == "region", ]
processed_He_results_regions$pop <- processed_He_results_regions$population


#write.csv(processed_He_results_regions, "He_results_regions_2022.csv", row.names = F)

#sample sizes for each population

sample_size_provinces <- combined_df_merged %>%
  group_by(year, province) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_provinces$pop <- sample_size_provinces$province
sample_size_provinces

combined_df_merged_nodry <- combined_df_merged %>%
  filter(!grepl("Dry", province))

sample_size_regions <- combined_df_merged_nodry %>%
  group_by(year, region) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_regions$pop <- sample_size_regions$region
sample_size_regions

# CALCULATE pairwise Fst 
# Fst = (Ht - Hs) / Ht [same as 1 - Hs / Ht]*;
# Ht is (n1*Hs1 + n2*Hs2) / n1+n2 [n is individuals of each pop]; 
# Hs = Hs1 + Hs2 / 2 [Hs1 and Hs2 are post_stat_mean (heterozygosity for each locus) calculated by moire for each pop]. calculated genome-wide average, just avergaae Hs1 and Hs2 of all loci
# THIS IS DONE FOR EACH LOCUS

# 1.- calculate Hs and Ht for each locus
# 2.- calculate Fst [(Ht - Hs) / Ht] through a Linear Mixed Model


### 1) FOR REGIONS

fts_input_regions <- merge(processed_He_results_regions[c("locus", "post_stat_mean", "pop")], sample_size_regions[c("unique_NIDA2_count", "pop")], by = c("pop"))

# Generate all possible combinations of unique values
unique_pops <- unique(fts_input_regions$pop)
combinations <- expand.grid(pop1 = unique_pops, pop2 = unique_pops)

# Create an empty dataframe to store the results
fst_results_df <- data.frame(matrix(ncol = 7, nrow = 1))
colnames(fst_results_df) <- c("pop1", "pop2", "Hs1", "Hs2", "n1", "n2", "locus")

# Loop through each locus
for (locus in unique(fts_input_regions$locus)) {
  locus_subset <- fts_input_regions[fts_input_regions$locus == locus,]
  
  # Loop through each combination of populations
  for (i in 1:nrow(combinations)) {
    pop1 <- as.character(combinations$pop1[i])
    pop2 <- as.character(combinations$pop2[i])
    
    Hs1 <- locus_subset[locus_subset$pop == pop1, ]$post_stat_mean
    n1 <- locus_subset[locus_subset$pop == pop1, ]$unique_NIDA2_count
    
    Hs2 <- locus_subset[locus_subset$pop == pop2, ]$post_stat_mean
    n2 <- locus_subset[locus_subset$pop == pop2, ]$unique_NIDA2_count
    
    # Check if both populations have data for this locus
    if (length(Hs1) > 0 && length(Hs2) > 0) {
      row <- c(pop1, pop2, Hs1, Hs2, n1, n2, locus)
    } else {
      # If one population does not have data, assign NA to the corresponding entry
      row <- c(pop1, pop2, NA, NA, NA, NA, locus)
    }
    
    fst_results_df <- rbind(fst_results_df, row)
  }
}

#formating
fst_results_df <- fst_results_df[-1,]
fst_results_df$Hs1 <- round(as.numeric(fst_results_df$Hs1), 3)
fst_results_df$Hs2 <- round(as.numeric(fst_results_df$Hs2),3)
fst_results_df$n1 <- round(as.numeric(fst_results_df$n1),3)
fst_results_df$n2 <- round(as.numeric(fst_results_df$n2),3)

if (nrow(fst_results_df) == nrow(combinations)*length(unique(fts_input_regions$locus))){
  print("Got all expected combinations for per locus pairwise Fst calculations")
}else{
  print("grab a coffee.")
}

fst_results_df <- fst_results_df[complete.cases(fst_results_df), ] # remove NA rows (those on which the amplicon wasn't present in both populations, hence Fst can't be calculated)

#calculate Hs for populations PER LOCUS
fst_results_df$Hs <- (fst_results_df$Hs1 + fst_results_df$Hs2) / 2
fst_results_df$Hs <- round(as.numeric(fst_results_df$Hs),3)

#Calculate Ht (uses sample sizes for each pop)
fst_results_df$Ht <- ((fst_results_df$n1 * fst_results_df$Hs1)+ (fst_results_df$n2 * fst_results_df$Hs2))/(fst_results_df$n1 + fst_results_df$n2)
fst_results_df$Ht <- round(as.numeric(fst_results_df$Ht),3)

# # Define a function to calculate FST
# calculate_FST <- function(data, indices) {
#   sampled_data <- data[indices, ]
#   
#   # Extract pre-calculated Hs and Ht values
#   Hs <- sampled_data$Hs
#   Ht <- sampled_data$Ht
#   
#   # Calculate FST
#   FST <- ((Ht - Hs) / Ht)
#   
#   return(FST)
# }
# 
# #for each locus, just in case it's needed
# fst_results_df$Fst <- (fst_results_df$Ht - fst_results_df$Hs)/fst_results_df$Ht
# 
# ### BOOTSTRAP: NON-PARAMETRIC APPROACH TO LME
# 
# set.seed(123) # For reproducibility
# n_permutations <- 1000
# observed_means <- aggregate(fst ~ comparison, data = FST_LLM, mean)
# 
# perm_means <- replicate(n_permutations, {
#   perm_data <- FST_LLM
#   perm_data$comparison <- sample(perm_data$comparison)
#   perm_agg <- aggregate(fst ~ comparison, data = perm_data, mean)
#   perm_agg$fst
# })
# 
# p_values <- sapply(1:nrow(observed_means), function(i) {
#   mean(observed_means$fst[i] <= perm_means[i, ])
# })
# 
# results <- data.frame(
#   comparison = observed_means$comparison,
#   observed_mean_fst = observed_means$fst,
#   p_value = p_values
# )
# 
# results



#llm (interchangeable with boostrat analysis)
library(nlme)
library(lmerTest)

FST_LLM <- as.data.frame(cbind(pop1 =fst_results_df$pop1,
                               pop2 = fst_results_df$pop2,
                               comparison = paste0(fst_results_df$pop1, "_",fst_results_df$pop2),
                               locus = fst_results_df$locus,
                               fst = as.numeric(fst_results_df$Fst)))


FST_LLM$fst <- as.numeric(FST_LLM$fst)

fst.model.region <- lme(fst ~ comparison,
                        random = ~ 1 | locus,
                        data = FST_LLM,
                        na.action = na.omit)

summary_data <- summary(fst.model.region)
summary_table <- as.data.frame(summary_data$tTable)
summary_table <- unique(summary_table)
summary_table$comparison <-  rownames(summary_table)

#dim(summary_table)

cis <- intervals(fst.model.region, which = "fixed")
cis <- as.data.frame(cis$fixed)
cis <- unique(cis)
cis$comparison <-  rownames(cis)

#dim(cis)

library(tidyr)

#merge estimates with CIs
final_table<- merge(cis, summary_table, by = c("comparison"))
final_table$comparison <- gsub("comparison", "", final_table$comparison)
final_table <- unique(merge(final_table, FST_LLM[c("pop1", "pop2", "comparison")], by = c("comparison")))

#heatmap
final_table$Fst_estimate<- round(final_table$est., 3)
final_table <- complete(final_table, pop1, pop2, fill = list(Fst_estimate = 0, `p-value` = 1))
final_table$Fst_estimate<- ifelse(final_table$Fst_estimate < 0, 0, final_table$Fst_estimate) #TURN NEGATIVE VALUES TO 0
final_table$label <- ifelse(final_table$`p-value` < 0.05 & final_table$Fst_estimate > 0 , paste0(final_table$Fst_estimate, "*"), as.character(final_table$Fst_estimate))

regions <- c("North", "Centre", "South") #ordered from north to south
#regions <- rev(regions)

final_table$pop1 <- factor(final_table$pop1, levels = regions)
final_table$pop2 <- factor(final_table$pop2, levels = regions)

heatmap_regions <- ggplot(final_table, aes(x = pop2, y = pop1, fill = Fst_estimate, label = label)) +
  geom_tile() +
  geom_text(color = "black") +
  scale_fill_gradient(low = "lightblue1", high = "orange", limits = c(min(final_table$Fst_estimate), max(final_table$Fst_estimate))) +  # Adjust scale limits
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "", y = "")

#significance
final_table$significance <- ifelse(final_table$`p-value` < 0.05, "p < 0.05", "not_signiff.")

final_table <- final_table %>%
  arrange(`est.`)

#keep fst > 0
final_table <- final_table[final_table$est. > 0.000001,]

# Create logical index to keep every other row
keep_rows <- seq(nrow(final_table)) %% 2 == 1

# Subset final_table to keep every other row
final_table <- final_table[keep_rows, ]

final_table <- final_table %>%
  mutate(comparison = factor(comparison, levels = comparison[order(est.)]))

final_table <- final_table[-nrow(final_table),]

anovap <- anova(fst.model.region, type = "marginal")
anovap_for_plot <- paste0("Anova's p = ", round(anovap$`p-value`[2], 3))

fst_regions <- ggplot(na.omit(final_table), aes(x = comparison, y = est., color = significance)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(title = "",
       x = "Pairwise comparisons",
       y = "Fst Estimate") +
  scale_color_manual(values = c("black", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_minimal()+
  annotate("text", x = Inf, y = -Inf, label = anovap_for_plot, hjust = 1.1, vjust = -1.1, size = 3, color = "black")

#fst_regions

#kw <- kruskal.test(fst ~ comparison, data = fst.model.region$data)
ggsave(paste0("fst_heatmap_regions_", SAMPLING, ".png"), heatmap_regions, width = 8, height = 6, bg = "white")
ggsave(paste0("fst_CI_regions_", SAMPLING, ".png"), fst_regions, width = 8, height = 6, bg = "white")
write.csv(final_table, paste0("Fst_regions_", SAMPLING, ".csv"), row.names = F)

### 2) FOR PROVINCES

fts_input_provinces <- merge(processed_He_results_provinces[c("locus", "post_stat_mean", "pop")], sample_size_provinces[c("unique_NIDA2_count", "pop")], by = c("pop"))

# Generate all possible combinations of unique values
unique_pops <- unique(fts_input_provinces$pop)
combinations <- expand.grid(pop1 = unique_pops, pop2 = unique_pops)

# Create an empty dataframe to store the results
fst_results_df <- data.frame(matrix(ncol = 7, nrow = 1))
colnames(fst_results_df) <- c("pop1", "pop2", "Hs1", "Hs2", "n1", "n2", "locus")

# Loop through each locus
for (locus in unique(fts_input_provinces$locus)) {
  locus_subset <- fts_input_provinces[fts_input_provinces$locus == locus,]
  
  # Loop through each combination of populations
  for (i in 1:nrow(combinations)) {
    pop1 <- as.character(combinations$pop1[i])
    pop2 <- as.character(combinations$pop2[i])
    
    Hs1 <- locus_subset[locus_subset$pop == pop1, ]$post_stat_mean
    n1 <- locus_subset[locus_subset$pop == pop1, ]$unique_NIDA2_count
    
    Hs2 <- locus_subset[locus_subset$pop == pop2, ]$post_stat_mean
    n2 <- locus_subset[locus_subset$pop == pop2, ]$unique_NIDA2_count
    
    # Check if both populations have data for this locus
    if (length(Hs1) > 0 && length(Hs2) > 0) {
      row <- c(pop1, pop2, Hs1, Hs2, n1, n2, locus)
    } else {
      # If one population does not have data, assign NA to the corresponding entry
      row <- c(pop1, pop2, NA, NA, NA, NA, locus)
    }
    
    fst_results_df <- rbind(fst_results_df, row)
  }
}

#formating
fst_results_df <- fst_results_df[-1,]
fst_results_df$Hs1 <- round(as.numeric(fst_results_df$Hs1), 3)
fst_results_df$Hs2 <- round(as.numeric(fst_results_df$Hs2),3)
fst_results_df$n1 <- round(as.numeric(fst_results_df$n1),3)
fst_results_df$n2 <- round(as.numeric(fst_results_df$n2),3)

if (nrow(fst_results_df) == nrow(combinations)*length(unique(fts_input_provinces$locus))){
  print("Got all expected combinations for per locus pairwise Fst calculations")
}else{
  print("grab a coffee.")
}


fst_results_df <- fst_results_df[complete.cases(fst_results_df), ] # remove NA rows (those on which the amplicon wasn't present in both populations, hence Fst can't be calculated)

#calculate Hs for populations PER LOCUS
fst_results_df$Hs <- (fst_results_df$Hs1 + fst_results_df$Hs2) / 2
fst_results_df$Hs <- round(as.numeric(fst_results_df$Hs),3)

#Calculate Ht (uses sample sizes for each pop)
fst_results_df$Ht <- ((fst_results_df$n1 * fst_results_df$Hs1)+ (fst_results_df$n2 * fst_results_df$Hs2))/(fst_results_df$n1 + fst_results_df$n2)
fst_results_df$Ht <- round(as.numeric(fst_results_df$Ht),3)

# Define a function to calculate FST
calculate_FST <- function(data, indices) {
  sampled_data <- data[indices, ]
  
  # Extract pre-calculated Hs and Ht values
  Hs <- sampled_data$Hs
  Ht <- sampled_data$Ht
  
  # Calculate FST
  FST <- ((Ht - Hs) / Ht)
  
  return(FST)
}

#for each locus, just in case it's needed
fst_results_df$Fst <- (fst_results_df$Ht - fst_results_df$Hs)/fst_results_df$Ht


#llm (interchangeable with boostrat analysis)
library(nlme)

FST_LLM <- as.data.frame(cbind(pop1 =fst_results_df$pop1,
                               pop2 = fst_results_df$pop2,
                               comparison = paste0(fst_results_df$pop1, "_",fst_results_df$pop2),
                               locus = fst_results_df$locus,
                               fst = as.numeric(fst_results_df$Fst)))

FST_LLM$fst <- as.numeric(FST_LLM$fst)

fst.model.region <- lme(fst ~ comparison,
                        random = ~ 1 | locus,
                        data = FST_LLM,
                        na.action = na.omit)

summary_data <- summary(fst.model.region)
summary_table <- as.data.frame(summary_data$tTable)
summary_table <- unique(summary_table)
summary_table$comparison <-  rownames(summary_table)

dim(summary_table)

cis <- intervals(fst.model.region, which = "fixed")
cis <- as.data.frame(cis$fixed)
cis <- unique(cis)
cis$comparison <-  rownames(cis)

dim(cis)

#merge estimates with CIs
final_table<- merge(cis, summary_table, by = c("comparison"))
final_table$comparison <- gsub("comparison", "", final_table$comparison)
final_table <- unique(merge(final_table, FST_LLM[c("pop1", "pop2", "comparison")], by = c("comparison")))

#heatmap
final_table$Fst_estimate<- round(final_table$est., 3)
final_table <- complete(final_table, pop1, pop2, fill = list(Fst_estimate = 0, `p-value` = 1))
final_table$Fst_estimate<- ifelse(final_table$Fst_estimate < 0, 0, final_table$Fst_estimate) #TURN NEGATIVE VALUES TO 0
final_table$label <- ifelse(final_table$`p-value` < 0.05 & final_table$Fst_estimate > 0 , paste0(final_table$Fst_estimate, "*"), as.character(final_table$Fst_estimate))


final_table$pop1 <- factor(final_table$pop1, levels = provinces)
final_table$pop2 <- factor(final_table$pop2, levels = provinces)

heatmap_provinces <- ggplot(final_table, aes(x = pop2, y = pop1, fill = Fst_estimate, label = label)) +
  geom_tile() +
  geom_text(color = "black") +
  scale_fill_gradient(low = "lightblue1", high = "orange", limits = c(min(final_table$Fst_estimate), max(final_table$Fst_estimate))) +  # Adjust scale limits
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "", y = "")

#significance
final_table$significance <- ifelse(final_table$`p-value` < 0.05, "p < 0.05", "not_signiff.")

final_table <- final_table %>%
  arrange(`est.`)

#keep fst > 0 
final_table <- final_table[final_table$est. > 0.000001,]

# Create logical index to keep every other row
keep_rows <- seq(nrow(final_table)) %% 2 == 1

# Subset final_table to keep every other row
final_table <- final_table[keep_rows, ]

final_table <- final_table %>%
  mutate(comparison = factor(comparison, levels = comparison[order(est.)]))

final_table <- final_table[-nrow(final_table),]

anovap <- anova(fst.model.region, type = "marginal")
anovap_for_plot <- paste0("Anova's p = ", round(anovap$`p-value`[2], 3))

fst_provinces <- ggplot(na.omit(final_table), aes(x = comparison, y = est., color = significance)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(title = "",
       x = "Pairwise comparisons",
       y = "Fst Estimate") +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  annotate("text", x = Inf, y = -Inf, label = anovap_for_plot, hjust = 1.1, vjust = -1.1, size = 3, color = "black")

#fst_regions

#kw <- kruskal.test(fst ~ comparison, data = fst.model.region$data)
ggsave(paste0("fst_heatmap_provinces_", SAMPLING, ".png"), heatmap_provinces, width = 8, height = 6, bg = "white")
ggsave(paste0("fst_CI_provinces_", SAMPLING, ".png"), fst_provinces, width = 8, height = 6, bg = "white")
write.csv(final_table, paste0("Fst_provinces_", SAMPLING, ".csv"), row.names = F)

