

library(ggplot2)
library(dplyr)
library(ggsignif)
library(tidyr)


regions <- c("South", "Centre", "North")

######################################################################################################
##### COMPARE eCOI #####

coi_2021 <- read.csv("coi_results_2021.csv")
coi_2022 <- read.csv("coi_results_2022.csv")

coi_all <- rbind(coi_2021, coi_2022)

coi_all$region <- factor(coi_all$region, levels = regions)

# Define a function to perform pairwise comparisons for each province
pairwise_comparisons <- function(data) {
  comparisons <- combn(unique(data$year), 2, simplify = FALSE)
  test_results <- lapply(comparisons, function(pair) {
    test <- wilcox.test(post_effective_coi_med ~ year, data = data[data$year %in% pair,])
    p_value <- test$p.value
    data.frame(region = unique(data$region), group1 = pair[1], group2 = pair[2], p_value = p_value)
  })
  do.call(rbind, test_results)
}

# Calculate significance values
significance_data <- coi_all %>%
  group_by(region) %>%
  do(pairwise_comparisons(.)) %>%
  ungroup()

# Prepare the data for geom_signif
significance_data <- significance_data %>%
  mutate(y_position = max(coi_all$post_effective_coi_med) + 0.5)  # Adjust y_position as needed

significance_data

# Plot with significance annotations
ECOI_PLOT <- ggplot(coi_all, aes(x = region, y = post_effective_coi_med, color = as.factor(year), fill = as.factor(year))) +
  geom_violin(width = 0.8, alpha = 0.4, position = position_dodge(width = 0.9)) +  # Use position_dodge to align
  geom_boxplot(width = 0.1, alpha = 0.7, position = position_dodge(width = 0.9), fill = "white") +  # Align boxplot with violin
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_discrete(name = "Region") +
  labs(x = "", y = "eMOI") +
  guides(color = FALSE, fill = FALSE)+
  scale_fill_manual(values = c("steelblue", "orange2"))+
  scale_color_manual(values = c("steelblue", "orange2"))

ECOI_PLOT

######################################################################################################
##### COMPARE %polyclonal #####

perc_polyclonal <-coi_all %>%
  group_by(region, year) %>%
  summarize(perc_polyclonal =  sum(post_effective_coi_med > 1.1) / length(post_effective_coi_med))

perc_polyclonal$region <- factor(perc_polyclonal$region, levels = regions)

#####33
library(dplyr)
library(purrr)
library(tidyr)

# Preparing data
coi_results_region_year <- coi_all %>%
  mutate(polyclonal = ifelse(post_effective_coi_med > 1.1, 1, 0))

# Function to perform pairwise Fisher's Exact Test for years within each region
pairwise_fisher_test_year <- function(data, year1, year2) {
  data1 <- data %>% filter(year == year1)
  data2 <- data %>% filter(year == year2)
  
  table <- matrix(c(
    sum(data1$polyclonal), nrow(data1) - sum(data1$polyclonal),
    sum(data2$polyclonal), nrow(data2) - sum(data2$polyclonal)
  ), nrow = 2, byrow = TRUE)
  
  fisher.test(table)$p.value
}

# Unique combinations of regions and years
regions <- unique(as.character(coi_results_region_year$region))

# Create an empty tibble to store the results
pairwise_results_year <- tibble()

# Loop through each region to perform pairwise comparisons between years
for (region in regions) {
  # Filter data for the current region
  data_region <- coi_results_region_year %>% filter(region == region)
  
  # Get the unique years for the current region
  years <- unique(data_region$year)
  
  # Generate all pairs of years for the current region
  year_pairs <- combn(years, 2, simplify = FALSE)
  
  # Perform pairwise Fisher's Exact Test for each pair of years
  results <- tibble(
    region = region,
    year1 = sapply(year_pairs, `[[`, 1),
    year2 = sapply(year_pairs, `[[`, 2),
    p.value = map2_dbl(sapply(year_pairs, `[[`, 1), sapply(year_pairs, `[[`, 2), 
                       ~pairwise_fisher_test_year(data_region, .x, .y))
  )
  
  # Append the results for the current region to the overall results
  pairwise_results_year <- bind_rows(pairwise_results_year, results)
}

# Adjust for multiple comparisons (e.g., Bonferroni correction)
pairwise_results_year <- pairwise_results_year %>%
  mutate(p.adjusted = p.adjust(p.value, method = "bonferroni"))

pairwise_results_year

########33


library(binom)

# Calculate proportion of polyclonal infections and confidence intervals
summary_data <- coi_all %>%
  mutate(polyclonal = ifelse(post_effective_coi_med > 1.1, 1, 0)) %>%
  group_by(year, region) %>%
  summarise(
    n = n(),
    n_poly = sum(polyclonal),
    prop_poly = mean(polyclonal)
  ) %>%
  mutate(
    ci_lower = binom.confint(n_poly, n, method = "wilson")$lower,
    ci_upper = binom.confint(n_poly, n, method = "wilson")$upper
  )

PERC_POLYCLONAL_PLOT <- ggplot(summary_data, aes(x = region, y = prop_poly, fill = as.factor(year))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.9)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_discrete(name = "Region") +
  labs(x = "", y = "%Polyclonal Samples") +
  guides(color = FALSE, fill = FALSE)+
  scale_fill_manual(values = c("steelblue", "orange2"))+
  scale_color_manual(values = c("steelblue", "orange2"))

PERC_POLYCLONAL_PLOT

######################################################################################################
##### COMPARE He #####

he_2021 <- read.csv("He_results_regions_2021.csv")
he_2021$year  <- 2021 
he_2022 <- read.csv("He_results_regions_2022.csv")
he_2022$year  <- 2022 

he_all <- rbind(he_2021, he_2022)

he_all<- he_all[he_all$geo == "region",]

he_all$population <- factor(he_all$population, levels = regions)


##################33
# library(emmeans)
# 
# # Fit separate linear models for each region
# region_models <- list()
# for (population in unique(he_all$population)) {
#   region_data <- he_all[he_all$population == population, ]
#   region_models[[population]] <- lm(post_stat_mean ~ year, data = region_data)
# }
# 
# # Perform pairwise comparisons between years within each region
# pairwise_comparisons <- list()
# for (population in unique(he_all$population)) {
#   region_model <- region_models[[population]]
#   pairwise_comparisons[[population]] <- emmeans(region_model, pairwise ~ year, adjust = "tukey")
# }
# 
# pairwise_comparisons
# ##################33

##################33
library(nlme)
library(emmeans)
library(lmerTest)

# Fit separate linear models for each region
region_models <- list()
for (population in unique(he_all$population)) {
  region_data <- he_all[he_all$population == population, ]
  region_models[[population]] <- lme(post_stat_mean ~ year,
                                     random = ~ 1 | locus, #locus is random effect
                                     data = region_data)
}


# Loop through each element in region_models
model_results <- list()

for (population in names(region_models)) {

  model <- region_models[[population]]

  summ <- summary(model)
  an <- anova(model)
  ci_mod <- intervals(model, which = "fixed")

  model_results[[paste0(population, "_summary")]] <- summ
  model_results[[paste0(population, "_anova")]] <- an
  model_results[[paste0(population, "_CI")]] <- ci_mod
}


# Perform pairwise comparisons between years within each region
pairwise_comparisons <- list()
for (population in unique(he_all$population)) {
  region_model <- region_models[[population]]
  pairwise_comparisons[[population]] <- emmeans(region_model, pairwise ~ year, adjust = "tukey")
}

pairwise_comparisons
# ##################33
# 
# #LMM: SALE SIGNIFICATIVO PARA TODAS LAS REGIONES
# 
# library(lme4)
# 
# # List to store results
# results_list <- list()
# 
# # Get unique regions
# unique_regions <- unique(he_all$population)
# 
# # Loop through each region
# for (region in unique_regions) {
#   # Subset data for the current region
#   region_data <- subset(he_all, region == region)
# 
#   # Fit the mixed effects model
#   mixed_model <- lmer(post_stat_med ~ year + (1 | locus), data = region_data)
# 
#   # Store the summary in the results list
#   results_list[[region]] <- summary(mixed_model)
# }
# 
# # Print the results
# for (region in unique_regions) {
#   cat("\nResults for region:", region, "\n")
#   print(results_list[[region]])
# }
# 

# ######## BOOTSTRAP: SALE SIGNIFICATIVO PARA TODAS LAS REGIONES
# 
# library(lme4) 
# 
# set.seed(123) # For reproducibility
# 
# # Number of bootstrap samples
# n_bootstrap <- 1000
# 
# # Get unique regions
# unique_regions <- unique(he_all$population)
# 
# # List to store bootstrap test results
# bootstrap_results_list <- list()
# 
# # Loop through each region
# for (region in unique_regions) {
#   # Subset data for the current region
#   region_data <- he_all[he_all$population == region,]
#   
#   # Bootstrap test
#   bootstrap_effects <- replicate(n_bootstrap, {
#     bootstrap_data <- region_data[sample(nrow(region_data), replace = TRUE), ]
#     bootstrap_model <- lmer(post_stat_med ~ year + (1 | locus), data = bootstrap_data)
#     fixef(bootstrap_model)["year"]
#   })
#   
#   ci_lower <- quantile(bootstrap_effects, 0.025)
#   ci_upper <- quantile(bootstrap_effects, 0.975)
#   
#   # Calculate p-value
#   p_value <- mean(bootstrap_effects <= 0) * 2 # Two-tailed test
#   
#   # Store results
#   bootstrap_results_list[[region]] <- list(ci_lower = ci_lower, ci_upper = ci_upper, p_value = p_value)
# }
# 
# # Define a function to perform pairwise comparisons for each province
# pairwise_comparisons <- function(data) {
#   comparisons <- combn(unique(data$year), 2, simplify = FALSE)
#   test_results <- lapply(comparisons, function(pair) {
#     test <- wilcox.test(post_stat_med ~ year, data = data[data$year %in% pair,])
#     p_value <- test$p.value
#     data.frame(region = unique(data$population), group1 = pair[1], group2 = pair[2], p_value = p_value)
#   })
#   do.call(rbind, test_results)
# }
# 
# # Calculate significance values
# significance_data <- he_all %>%
#   group_by(population) %>%
#   do(pairwise_comparisons(.)) %>%
#   ungroup()
# 
# # Prepare the data for geom_signif
# significance_data <- significance_data %>%
#   mutate(y_position = max(he_all$post_stat_med) + 0.5)  # Adjust y_position as needed
# 
# significance_data


# Extract EMMs and CIs into a dataframe
emm_data <- lapply(names(pairwise_comparisons), function(region) {
  emm <- as.data.frame(pairwise_comparisons[[region]]$emmeans)
  emm$region <- region
  return(emm)
})

emm_data <- bind_rows(emm_data)

colnames(emm_data) <- c("year", "emmean", "SE", "df", "lower.CL", "upper.CL", "region")

emm_data$region <- factor(emm_data$region, levels = regions)

emm_p <- lapply(names(pairwise_comparisons), function(region) {
  emm <- as.data.frame(pairwise_comparisons[[region]]$contrasts[1])
  emm$region <- region
  return(emm)
})

# Plotting the data
HE_PLOT<- ggplot(emm_data, aes(x = factor(region), y = emmean, color = factor(year))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, position = position_dodge(width = 0.5),linewidth = 1) +
  geom_text(data = emm_p[[1]], aes(label = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", "")))),
            x = "Centre", y = max(emm_data$upper.CL)+0.005, vjust = -0.5, color = "black", size = 6) +
  geom_text(data = emm_p[[2]], aes(label = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", "")))),
            x = "North", y = max(emm_data$upper.CL)+0.005, vjust = -0.5, color = "black", size = 6) +
  geom_text(data = emm_p[[3]], aes(label = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", "")))),
            x = "South", y = max(emm_data$upper.CL)+0.005, vjust = -0.5, color = "black", size = 6) +
  labs(title = "",
       x = "",
       y = "He estimate",
       color = "Year") +
  theme_minimal() +
  guides(color = guide_legend(title = "Year"))+
  ylim(min(emm_data$emmean)-0.05, max(emm_data$emmean)+0.05)+
  guides(color = FALSE, fill = FALSE)+
  scale_fill_manual(values = c("steelblue", "orange2"))+
  scale_color_manual(values = c("steelblue", "orange2"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


HE_PLOT

# ######################################################################################################
# ### COMPARE CONNECTIVITY ###
# 
# conn_2021 <- read.csv("connectivity_2021.csv")
# conn_2022 <- read.csv("connectivity_2022.csv")
# 
# conn_2021$year <- 2021
# conn_2022$year <- 2022
# 
# conn_all <- rbind(conn_2021, conn_2022)
# 
# # Define a function to perform pairwise comparisons for each province
# pairwise_comparisons <- function(data) {
#   comparisons <- combn(unique(data$year), 2, simplify = FALSE)
#   test_results <- lapply(comparisons, function(pair) {
#     test <- wilcox.test(estimate ~ year, data = data[data$year %in% pair,])
#     p_value <- test$p.value
#     data.frame(region = unique(data$conn_regions), group1 = pair[1], group2 = pair[2], p_value = p_value)
#   })
#   do.call(rbind, test_results)
# }
# 
# # Calculate significance values
# significance_data <- conn_all %>%
#   group_by(conn_regions) %>%
#   do(pairwise_comparisons(.)) %>%
#   ungroup()
# 
# # Prepare the data for geom_signif
# significance_data <- significance_data %>%
#   mutate(y_position = max(conn_all$estimate) + 0.5)  # Adjust y_position as needed
# 
# significance_data
# 
# 
# # Plot with significance annotations
# CONN_PLOT <- ggplot(conn_all, aes(x = conn_regions, y = estimate, color = as.factor(year)))+
#   geom_boxplot() +
#   theme_minimal() +
#   scale_x_discrete(labels = function(x) gsub("\\..*", "", x)) +  # Remove year from x-axis labels
#   labs(x = "Region", y = "IBD", color = "Year", fill = "Year", title = "Connectivity") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))



######################################################################################################
### COMPARE 1-FWS ###

fws_2021 <- read.csv("He_and_Fws_2021.csv")
fws_2022 <- read.csv("He_and_Fws_2022.csv")

fws_all <- rbind(fws_2021, fws_2022)
fws_all <- fws_all[!is.na(fws_all$fws_region),]

fws_all$region <- factor(fws_all$region, levels = regions)

# remove 5% least variant amplicons to avoid issues when calculating Fws
amps_q95<- round(length(unique(fws_all$locus))*0.95)

he_amps<- fws_all %>%
  group_by(locus) %>%
  summarize(He_region = mean(He_region))%>%
  arrange(desc(He_region))

he_amps_topHe <- he_amps[1:amps_q95,]

fws_all<- fws_all[fws_all$locus %in% he_amps_topHe$locus,]

fws_all <- fws_all[,c("NIDA2", "locus", "year", "region", "fws_region")] 

fws_all <- distinct(fws_all)

fws_all <- fws_all %>%
  group_by(NIDA2, year, region) %>%
  summarize(fws_region = mean(fws_region))


# Define a function to perform pairwise comparisons for each province
pairwise_comparisons <- function(data) {
  comparisons <- combn(unique(data$year), 2, simplify = FALSE)
  test_results <- lapply(comparisons, function(pair) {
    test <- wilcox.test(fws_region ~ year, data = data[data$year %in% pair,])
    p_value <- test$p.value
    data.frame(region = unique(data$region), group1 = pair[1], group2 = pair[2], p_value = p_value)
  })
  do.call(rbind, test_results)
}

# Calculate significance values
significance_data <- fws_all %>%
  group_by(region) %>%
  do(pairwise_comparisons(.)) %>%
  ungroup()

# Prepare the data for geom_signif
significance_data <- significance_data %>%
  mutate(y_position = max(fws_all$fws_region) + 0.5)  # Adjust y_position as needed

significance_data


FWS_PLOT <- ggplot(fws_all, aes(x = region, y = fws_region, color = as.factor(year), fill = as.factor(year)))+
  geom_violin(width = 0.8, alpha = 0.4, position = position_dodge(width = 0.9)) +  # Use position_dodge to align
  geom_boxplot(width = 0.1, alpha = 0.7, position = position_dodge(width = 0.9), fill = "white") +  # Align boxplot with violin
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Genomewide 1-Fws") +
  guides(color = FALSE, fill = FALSE)+
  scale_fill_manual(values = c("steelblue", "orange2"))+
  scale_color_manual(values = c("steelblue", "orange2"))


FWS_PLOT

###################################################################################################
# FINAL FIGURES

library(gridExtra)

combined_plot <- grid.arrange(
  ECOI_PLOT, PERC_POLYCLONAL_PLOT, FWS_PLOT, HE_PLOT,
  ncol = 2, nrow = 2
)

combined_plot


ggsave("comparison_plots_22222.png", combined_plot, width = 10, height = 7)
