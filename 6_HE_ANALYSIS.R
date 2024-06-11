
library(ggsignif)
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
  
  provinces <- c("Niassa", "Nampula", "Zambezia", "Manica", "Inhambane", "Maputo") #ordered from north to south
  province_colors <- c(Niassa = "firebrick4", Nampula = "indianred1", Zambezia = "darkgreen", Manica = "green", Inhambane = "cornflowerblue", Maputo = "deepskyblue")
  
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
processed_He_results$population <- gsub(paste0("_",SAMPLING, "_MOIRE-RESULTS_FOR_ALLELE_FREQS"), "", processed_He_results$population)

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




## linear mixed models to assess difference in He (from Nanna's scripts)
# LLM (EACH SITE AS REFERENCE)
library(nlme)

# Define the population levels
regions <- c("North", "Centre", "South")


###
# Function to fit LLM and collect results for each reference level
fit_and_collect_results <- function(reference_population, pop = "region") {
  
  # Subset data and set reference population
  he_province <- processed_He_results[processed_He_results$geo == pop,] %>%
    mutate(population = factor(population, levels = population_levels)) %>%
    mutate(population = relevel(population, reference_population))
  
  # Fit LLM
  he.model.province <- lme(post_stat_med ~ population,
                           random = ~ 1 | locus,
                           data = he_province,
                           na.action = na.omit)
  
  # Normality check
  residuals <- resid(he.model.province)
  shapiro_result <- shapiro.test(residuals)  # Shapiro-Wilk test for normality
  
  # Extract and return summary statistics
  summary_data <- summary(he.model.province)
  aic <- AIC(logLik(he.model.province))
  t_table <- summary_data$tTable
  ci_mod <- intervals(he.model.province, which = "fixed")
  p_values <- anova(he.model.province, type = "marginal")$'p-value'[2]
  
  # Extract CI for reference population
  ci_reference <- ci_mod$fixed[1,]
  ci_reference <- t(as.data.frame(ci_reference))
  rownames(ci_reference) <- reference_population
  
  # Extract CIs for other populations
  ci_others <- t(sapply(2:length(population_levels), function(i) ci_reference + ci_mod$fixed[i,]))
  rownames(ci_others) <- population_levels[!population_levels %in% reference_population]
  
  ci_all <- rbind(ci_reference, ci_others)
  
  #merge tables
  t_table <- cbind(t_table, ci_all)
  rownames(t_table) <- rownames(ci_all)
  t_table <- as.data.frame(t_table)
  t_table$reference_pop <- reference_population
  t_table$compared_pop <- rownames(t_table)
  
  # Combine results
  results <- list(reference_population = reference_population,
                  t_table = t_table,
                  shapiro_result = shapiro_result)
  
  return(results)
}

###

analyze_results <- function(results_list, pop = "region") {
  # edit list
  names(results_list) <- population_levels
  
  results_list <- lapply(results_list, function(sublist) {
    # Remove the reference_population element from each sublist
    sublist$reference_population <- NULL
    return(sublist)
  })
  
  
  llm_tables_list <- list()
  
  # Extract t_table from each element in results_list
  for (population_name in names(results_list)) {
    t_table <- results_list[[population_name]]$t_table
    llm_tables_list[[population_name]] <- t_table
  }
  
  # Combine all t_tables into a single dataframe
  llm_tables_all <- do.call(rbind, llm_tables_list)
  
  llm_tables_all$lower <- round(llm_tables_all$lower, 3)
  llm_tables_all$est. <- round(llm_tables_all$est., 3)
  llm_tables_all$upper <- round(llm_tables_all$upper, 3)
  
  llm_tables_all <- llm_tables_all[llm_tables_all$reference_pop != llm_tables_all$compared_pop,]
  
  #plot He and 95% CI for each pop
  CIs <- llm_tables_all[c("compared_pop", "lower", "est.", "upper")]
  rownames(CIs) <- NULL
  
  CIs <- unique(CIs)
  
  CIs$compared_pop <- factor(CIs$compared_pop, levels = population_levels)
  
  #put zero to self comparisons
  #llm_tables_all[llm_tables_all$`p-value` < 0.0000001, ]$`p-value` <- 0
  
  #significant pairwise comparions
  llm_tables_significant_pariwise<- llm_tables_all[llm_tables_all$reference_pop != llm_tables_all$compared_pop &  llm_tables_all$`p-value`< 0.05,]
  llm_tables_significant_pariwise <- subset(llm_tables_significant_pariwise, !duplicated(t(apply(llm_tables_significant_pariwise[c("compared_pop", "reference_pop")], 1, sort))))
  
  alpha <- 0.05
  sig_levels <- c("***", "**", "*", "")
  
  # Define significance labels based on p-values
  llm_tables_significant_pariwise$significance <- cut(llm_tables_significant_pariwise$`p-value`, 
                                                      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                                      labels = sig_levels,
                                                      include.lowest = TRUE)
  
  #AIC AND ANOVA'S P-VALUE (it's the same for every comparison)
  
  # Subset data and set reference population
  he_province <- processed_He_results[processed_He_results$geo == pop ,] %>%
    mutate(population = factor(population, levels = population_levels))
  
  he.model.province <- lme(post_stat_med ~ population,
                           random = ~ 1 | locus,
                           data = he_province,
                           na.action = na.omit)
  
  anovap <- anova(he.model.province, type = "marginal")
  anovap_for_plot <- paste0("Anova's p = ", round(anovap$`p-value`[2], 3))
  
  aicval <- AIC(logLik(he.model.province))
  
  #kw <- kruskal.test(post_stat_med ~ population, data = he.model.province$data)
  
  plotci <- ggplot(CIs, aes(x = compared_pop, y = est.)) +
    geom_point(color = "black") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "black") +
    labs(title = "",
         x = "",
         y = "He Estimate") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    annotate("text", x = Inf, y = -Inf, label = anovap_for_plot, hjust = 1.1, vjust = -1.1, size = 3, color = "black")
  
  
  # Combine results
  results <- list(plotCI = plotci,
                  llm_tables_significant_pariwise = llm_tables_significant_pariwise,
                  anovap = anovap,
                  AIC = aicval) #,
                  #KW = kw)
  
  return(results)
}
###

# PROVINCES
population_levels <- provinces
results_list_provinces <- lapply(provinces, fit_and_collect_results, pop = "province")
results_provinces <- analyze_results(results_list_provinces, pop = "province")

##############################
#add significance to plot
plotp <- results_provinces[[1]]
sigs <- results_provinces$llm_tables_significant_pariwise
anovap_for_plot <- paste0("Anova's p = ", round(results_provinces$anovap$`p-value`[2], 3))


# Create a list of comparisons for ggsignif
comparisons <- lapply(1:nrow(sigs), function(i) {
  c(as.character(sigs$reference_pop[i]), as.character(sigs$compared_pop[i]))
})

final_plot <- ggplot(plotp$data, aes(x = compared_pop, y = est.)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  theme_minimal() +
  labs(title = "",
       x = "Province",
       y = "He estimate") +
  geom_signif(
    comparisons = comparisons,
    annotations = sigs$significance,
    y_position = seq(max(plotp$data$upper) + 0.01, by = 0.01, length.out = nrow(sigs)),
    tip_length = 0.01
  )+
  annotate("text", x = Inf, y = -Inf, label = anovap_for_plot, hjust = 1.1, vjust = -1.1, size = 3, color = "black")
################################


write.csv(results_provinces$llm_tables_significant_pariwise, paste0("lmm_tables_significant_pariwise_provinces_", SAMPLING, ".csv"), row.names = F)
ggsave(paste0("province_He_lmm_", SAMPLING, ".png"), final_plot, width = 8, height = 6, bg = "white")



#PREGIONS
population_levels <- regions
results_list_regions <- lapply(regions, fit_and_collect_results, pop = "region")
results_regions <- analyze_results(results_list_regions, pop = "region")

##############################
#add significance to plot
plotp <- results_regions[[1]]
sigs <- results_regions$llm_tables_significant_pariwise
anovap_for_plot <- paste0("Anova's p = ", round(results_regions$anovap$`p-value`[2], 3))

# Create a list of comparisons for ggsignif
comparisons <- lapply(1:nrow(sigs), function(i) {
  c(as.character(sigs$reference_pop[i]), as.character(sigs$compared_pop[i]))
})

final_plot <- ggplot(plotp$data, aes(x = compared_pop, y = est.)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  theme_minimal() +
  labs(title = "",
       x = "Region",
       y = "He estimate") +
  geom_signif(
    comparisons = comparisons,
    annotations = sigs$significance,
    y_position = seq(max(plotp$data$upper) + 0.01, by = 0.01, length.out = nrow(sigs)),
    tip_length = 0.01
  )
################################

write.csv(results_regions$llm_tables_significant_pariwise, paste0("lmm_tables_significant_pariwise_regions_", SAMPLING, ".csv"), row.names = F)
ggsave(paste0("region_He_lmm_", SAMPLING, ".png"), final_plot, width = 8, height = 6, bg = "white")
