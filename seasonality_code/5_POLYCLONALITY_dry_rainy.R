library(dplyr)
library(haven)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)

SAMPLING <- 2022 # 2021 or 2022

####################################################

#metadata
db <- read_dta('../DrugRes_2021_2022_DB_ALLDATA_1Jan2024.dta')

combined_df_merged <- readRDS(paste0("combined_df_merged_", SAMPLING, "_only.RDS")) 
combined_df_merged<- combined_df_merged[combined_df_merged$province %in% c("Maputo_Rainy", "Maputo_Dry", "Manica_Rainy", "Manica_Dry"),]

#split province and season
combined_df_merged<- combined_df_merged %>%
  separate(province, into = c("province", "season"), sep = "_")

#import moire data
mcmc_results <- readRDS(paste0("all_samples_complete_filtered_MOIRE-RESULTS_", SAMPLING, "_only_FOR_MOI_usar_DRY.RDS"))

provinces <- c("Manica","Maputo") #ordered from north to south
# province_colors <- c(Manica = "forestgreen", Maputo = "cornflowerblue")

#######################################################
# 7.- Present MOI/eMOI results overall and means per province and region for each year
#######################################################

eff_coi <- moire::summarize_effective_coi(mcmc_results)
naive_coi <- moire::summarize_coi(mcmc_results)

# Merge the summaries by sample_id
coi_results <- merge(eff_coi, naive_coi, by = "sample_id")

# label mono and poly infections. NOTE: "proportion of polyclonal infections (eMOI>1.1)" from Nanna's manuscript
coi_results$polyclonal_from_ecoi_med <- ifelse(coi_results$post_effective_coi_med > 1.1, "polyclonal", "monoclonal")

#merge with categorical variables
colnames(coi_results)[1] <- "NIDA2"
coi_results <- merge(coi_results, db[c("NIDA2", "province", "region", "seasonality")], by="NIDA2")

coi_results$province <- gsub(" ", "_", coi_results$province) # for cabo delgado

# MANAGE SEASONALITY
seasonality_factor <- as_factor(coi_results$seasonality)
seasonality_recode <- recode(seasonality_factor, "1" = "Rainy", "2" = "Dry")
seasonality_char <- as.character(seasonality_recode)
coi_results$seasonality <- seasonality_char

coi_results$province <- ifelse(
  is.na(coi_results$seasonality),
  coi_results$province,
  paste0(coi_results$province, "_", coi_results$seasonality)
)

# coi_results_region$year <- 2022
# write.csv(coi_results_region, "coi_results_2022.csv",  row.names = F)


coi_results<- coi_results[coi_results$province %in% c("Maputo_Rainy", "Maputo_Dry", "Manica_Rainy", "Manica_Dry"),]

coi_results <- coi_results %>%
  separate(province, into = c("province", "season"), sep = "_")


polyclonal_percentage_province <- coi_results %>%
  group_by(province, season) %>%
  summarise(polyclonal_percentage_province = mean(polyclonal_from_ecoi_med == "polyclonal") * 100) %>%
  ungroup()

regions <- c("North", "Centre", "South")

coi_results$province <- factor(coi_results$province, levels = provinces)
polyclonal_percentage_province$province <- factor(polyclonal_percentage_province$province, levels = provinces)


# Convert the 'season' column to a factor for better visualization
coi_results$season <- as.factor(coi_results$season)

# Retrieve the default ggplot2 color palette
default_colors <- scales::hue_pal()(2)

# Invert the order of the colors
inverted_colors <- default_colors

# Define the colors for the seasons, matching the order of levels
season_colors <- setNames(inverted_colors, c("Dry", "Rainy"))

coi_results$season <- factor(coi_results$season, levels = c("Rainy", "Dry"))
coi_results$province <- factor(coi_results$province, levels = c("Maputo", "Manica"))

ecoi_plot <- ggplot(coi_results, aes(x = province, y = post_effective_coi_med, fill = season)) +
  geom_violin(width = 0.5, aes(color = season), alpha = 0.4) +
  geom_boxplot(width = 0.5, aes(color = season), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  labs(title = "",
       x = "",
       y = "eCOI")+
  scale_fill_manual(values = season_colors)+
  scale_color_manual(values = season_colors)

# Function to run Kruskal-Wallis test for each province
kruskal_results <- coi_results %>%
  group_by(province) %>%
  summarise(kruskal_p_value = kruskal.test(post_effective_coi_med ~ season)$p.value)

# Print the results
print(kruskal_results)

ggsave(paste0("ecoi_dry_rainy_violin_", SAMPLING,".png"), ecoi_plot, width = 8, height = 6, bg = "white")






## %POLYCLONAL

# f <- ggplot(polyclonal_percentage_province, aes(x = province, y = polyclonal_percentage_province,  fill = season))+
#   geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
#   labs(x = "", y = "%Polyclonal Infections") +
#   #facet_wrap(~province, scales = "fixed", nrow = 1) +
#   #scale_fill_manual(values = c("2022" = "orange")) + 
#   theme_minimal()+
#   #scale_fill_manual(values = province_colors)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggsave(paste0("perc_polyclonal_dry_rainy_", SAMPLING, ".png"), f, width = 8, height = 6, bg = "white")



library(tidyr)
library(purrr)
library(binom)
library(ggsignif)


# Preparing data
coi_results <- coi_results %>%
  mutate(polyclonal = ifelse(post_effective_coi_med > 1.1, 1, 0))


# Unique provinces
provinces <- unique(coi_results$province)

# Initialize an empty data frame to store results
fisher_results <- data.frame(province = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop over each province
for (province in provinces) {
  
  province_data <- coi_results %>% filter(province == !!province)

  contingency_table <- table(province_data$season, province_data$polyclonal)
  
  # Perform Fisher's Exact Test
  fisher_test_result <- fisher.test(contingency_table)

  fisher_results <- rbind(fisher_results, data.frame(province = province, p_value = fisher_test_result$p.value))
}

fisher_results


# Calculate proportion of polyclonal infections and confidence intervals
summary_data <- coi_results %>%
  mutate(polyclonal = ifelse(post_effective_coi_med > 1.1, 1, 0)) %>%
  group_by(province, season) %>%
  summarise(
    n = n(),
    n_poly = sum(polyclonal),
    prop_poly = mean(polyclonal)
  ) %>%
  mutate(
    ci_lower = binom.confint(n_poly, n, method = "wilson")$lower,
    ci_upper = binom.confint(n_poly, n, method = "wilson")$upper
  )

summary_data
summary_data$season <- factor(summary_data$season, levels = c("Rainy", "Dry"))
summary_data$province <- factor(summary_data$province, levels = c("Maputo", "Manica"))

# Retrieve the default ggplot2 color palette
default_colors <- scales::hue_pal()(2)

# Invert the order of the colors
inverted_colors <- default_colors

# Define the colors for the seasons, matching the order of levels
season_colors <- setNames(inverted_colors, c("Dry", "Rainy"))

# Create the plot
e <- ggplot(summary_data, aes(x = province, y = prop_poly, fill = season)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = "",
       x = "",
       y = "Proportion of Polyclonal Infections",
       fill = "Season") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = season_colors)

e

ggsave(paste0("perc_polyclonal_dry_rainy_", SAMPLING, ".png"), e, width = 8, height = 6, bg = "white")


# 
# # Prepare significance labels
# pairwise_results <- pairwise_results %>%
#   mutate(label = case_when(
#     p.adjusted < 0.001 ~ "***",
#     p.adjusted < 0.01 ~ "**",
#     p.adjusted < 0.05 ~ "*",
#     TRUE ~ ""
#   ))
# pairwise_results <- pairwise_results[pairwise_results$p.adjusted < 0.05,]
# 
# # Add significance annotations only for significant results
# for (i in 1:nrow(pairwise_results)) {
#   if (pairwise_results$label[i] != "") {
#     e <- e + geom_signif(
#       comparisons = list(c(pairwise_results$region1[i], pairwise_results$region2[i])),
#       annotations = pairwise_results$label[i],
#       map_signif_level = TRUE,
#       y_position = max(summary_data$prop_poly) + 0.05 * (i + 1)
#     )
#   }
# }


###
# coi results for Simone
coi_for_db <- coi_results[c("NIDA2", "naive_coi", "post_effective_coi_med")]

write.csv(coi_for_db, paste0("coi_for_db_dry_rainy", SAMPLING, ".csv"), row.names = F)

#plot(coi_for_db$naive_coi, coi_for_db$post_effective_coi_med)

scatter_plot <- ggplot(coi_for_db, aes(x = naive_coi, y = post_effective_coi_med)) +
  geom_point(alpha = 0.75) +
  labs(x = "Naive COI", y = "Post Effective COI Median") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("scatter_plot_dry_rainy",SAMPLING,".png"), scatter_plot, width = 8, height = 5, bg = "white")

