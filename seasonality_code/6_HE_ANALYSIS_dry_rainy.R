
library(ggsignif)
library(ggplot2)
library(dplyr)

SAMPLING <- 2022 # 2021 or 2022

####################################################

#metadata
db <- read_dta('../DrugRes_2021_2022_DB_ALLDATA_1Jan2024.dta')

combined_df_merged <- readRDS(paste0("combined_df_merged_", SAMPLING, "_only.RDS")) 
combined_df_merged<- combined_df_merged[combined_df_merged$province %in% c("Maputo_Rainy", "Maputo_Dry", "Manica_Rainy", "Manica_Dry"),]

#split province and season
combined_df_merged<- combined_df_merged %>%
  separate(province, into = c("province", "season"), sep = "_")

provinces <- c("Manica_Dry", "Manica_Rainy", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
# province_colors <- c(Manica = "forestgreen", Maputo = "cornflowerblue")


#######################################################
# 9.- He and Fws results 
#######################################################

# 1) calculate heterozygosity of the population (He); pop = province, region
#import everything into lists
rds_files <- list.files(path = ".", pattern = "\\MOIRE-RESULTS_FOR_ALLELE_FREQS.RDS$", full.names = TRUE)


rds_files <- rds_files[!rds_files %in% paste0("all_samples_complete_filtered_MOIRE-RESULTS_", SAMPLING, "_only_FOR_MOI_use_DRY.RDS")] 

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

processed_He_results<- processed_He_results %>%
  separate(population, into = c("province", "season"), sep = "_")


# #quick viz
# ggplot(processed_He_results, aes(x = province, y = post_stat_med, fill = season)) +
#   geom_violin(width = 0.5, aes(color = season), alpha = 0.4) +
#   geom_boxplot(width = 0.5, aes(color = season), fill = "white", alpha = 0.4) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 11), 
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )+
#   labs(title = "",
#        x = "",
#        y = "He Median")




library(nlme)

maputo <- processed_He_results[processed_He_results$province == "Maputo",]
manica <- processed_He_results[processed_He_results$province == "Manica",]

#maputo
he.model.maputo <- lme(post_stat_med ~ season,
                         random = ~ 1 | locus,
                         data = maputo,
                         na.action = na.omit)

summary_data <- summary(he.model.maputo)
aic <- AIC(logLik(he.model.maputo))
t_table <- summary_data$tTable
ci_mod <- intervals(he.model.maputo, which = "fixed")
p_values <- anova(he.model.maputo, type = "marginal")
anovas_p_maputo = p_values$`p-value`[2]

# Extract CI for reference population
ci_reference <- ci_mod$fixed[1,]
ci_reference <- t(as.data.frame(ci_reference))
rownames(ci_reference) <- "Maputo_Dry"

# Extract CIs for other populations
ci_others <- ci_reference + ci_mod$fixed[2,]
rownames(ci_others)<- "Maputo_Rainy"
ci__maputo <- rbind(ci_reference, ci_others)



#manica
he.model.manica <- lme(post_stat_med ~ season,
                       random = ~ 1 | locus,
                       data = manica,
                       na.action = na.omit)

summary_data <- summary(he.model.manica)
aic <- AIC(logLik(he.model.manica))
t_table <- summary_data$tTable
ci_mod <- intervals(he.model.manica, which = "fixed")
p_values <- anova(he.model.manica, type = "marginal")
anovas_p_manica = p_values$`p-value`[2]

# Extract CI for reference population
ci_reference <- ci_mod$fixed[1,]
ci_reference <- t(as.data.frame(ci_reference))
rownames(ci_reference) <- "Manica_Dry"

# Extract CIs for other populations
ci_others <- ci_reference + ci_mod$fixed[2,]
rownames(ci_others)<- "Manica_Rainy"
ci__manica <- rbind(ci_reference, ci_others)


#merge everything:
CI_ALL_all<- as.data.frame(rbind(ci__manica, ci__maputo))
CI_ALL_all$province <- rownames(CI_ALL_all)

CI_ALL_all<- CI_ALL_all %>%
  separate(province, into = c("province", "season"), sep = "_")

anovasp <- list(anovas_p_manica, anovas_p_maputo)

library(ggsignif)
#PLOT
dodge <- position_dodge(width = 0.5)


# Function to convert p-values to significance symbols
get_significance_symbol <- function(pval) {
  if (pval < 0.001) {
    return("***")
  } else if (pval < 0.01) {
    return("**")
  } else if (pval < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# Create the data frame for the text annotations with significance symbols
pval_annotations <- data.frame(
  province = c("Manica", "Maputo"),
  y = c(0.63, 0.63),
  pval = c(anovasp[[1]], anovasp[[2]])
)

# Add significance symbols based on p-values
pval_annotations$significance <- sapply(pval_annotations$pval, get_significance_symbol)
pval_annotations <- pval_annotations[rev(seq_len(nrow(pval_annotations))), ]

CI_ALL_all$season <- factor(CI_ALL_all$season, levels = c("Rainy", "Dry"))
CI_ALL_all$province <- factor(CI_ALL_all$province, levels = c("Maputo", "Manica"))

# Retrieve the default ggplot2 color palette
default_colors <- scales::hue_pal()(2)

# Invert the order of the colors
inverted_colors <- default_colors

# Define the colors for the seasons, matching the order of levels
season_colors <- setNames(inverted_colors, c("Dry", "Rainy"))

dry_rainy_plot <- ggplot(CI_ALL_all, aes(x = province, y = est., fill = season, color = season)) +
  geom_point(position = dodge) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = dodge) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "",
       x = "",
       y = "He estimate") +
  # Add text annotations for p-values with significance symbols
  geom_text(data = pval_annotations,
            aes(x = province, y = y, label = significance), 
            inherit.aes = FALSE,
            color = "black",
            size = 4) +
  # Add significance lines
  geom_segment(data = pval_annotations,
               aes(x = as.numeric(factor(province)) - 0.2, 
                   xend = as.numeric(factor(province)) + 0.2, 
                   y = y  - 0.003, 
                   yend = y - 0.003), 
               inherit.aes = FALSE, 
               color = "black")+
  scale_fill_manual(values = season_colors)+
  scale_color_manual(values = season_colors)

dry_rainy_plot

# # Function to fit the linear mixed model for each province
# fit_lmm_for_province <- function(data, province) {
#   model <- lmer(post_stat_med ~ season + (1 | locus), data = data)
#   summary(model)
# }
# 
# library(lmerTest)
# lmer_results <- processed_He_results %>%
#   group_by(province) %>%
#   do(model_summary = fit_lmm_for_province(., .$province))
# 
# lmer_results$model_summary[2]


ggsave(paste0("dry_rainy_He_lmm_", SAMPLING, ".png"), dry_rainy_plot, width = 8, height = 6, bg = "white")
