library(dplyr)
library(haven)
library(ggplot2)
library(reshape2)
library(ggpubr)

SAMPLING <- 2022 # 2021 or 2022

####################################################

#metadata
db <- read_dta('DrugRes_2021_2022_DB_ALLDATA_1Jan2024.dta')

combined_df_merged <- readRDS(paste0("combined_df_merged_", SAMPLING, "_only.RDS")) 

mcmc_results <- readRDS(paste0("all_samples_complete_filtered_MOIRE-RESULTS_", SAMPLING, "_only_FOR_MOI.RDS"))

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

# remove DRY season pops from region analysis
coi_results <- coi_results[!(coi_results$province %in% c("Maputo_Dry", "Manica_Dry")), ]
coi_results_region <- coi_results[!(coi_results$province %in% c("Maputo_Dry", "Manica_Dry")), ] 

# coi_results_region$year <- 2022
# write.csv(coi_results_region, "coi_results_2022.csv",  row.names = F)

# % polyclonal infections on each province and region
polyclonal_percentage_region <- coi_results_region %>%
  group_by(region) %>%
  summarise(polyclonal_percentage_region = mean(polyclonal_from_ecoi_med == "polyclonal") * 100) %>%
  ungroup()

polyclonal_percentage_province <- coi_results %>%
  group_by(province) %>%
  summarise(polyclonal_percentage_province = mean(polyclonal_from_ecoi_med == "polyclonal") * 100) %>%
  ungroup()

regions <- c("North", "Centre", "South")

coi_results$province <- factor(coi_results$province, levels = provinces)
coi_results$region <- factor(coi_results$region, levels = regions)
coi_results_region$region <- factor(coi_results_region$region, levels = regions)
polyclonal_percentage_region$region <- factor(polyclonal_percentage_region$region, levels = regions)
polyclonal_percentage_province$province <- factor(polyclonal_percentage_province$province, levels = provinces)


# #NAIVE COI
#
# a <- ggplot(coi_results, aes(x = naive_coi, fill = province)) +
#   geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
#   facet_wrap(~ province , scales = "fixed", nrow = 1) + 
#   labs(title = "",
#        x = "Naive COI",
#        y = "Frequency",
#        fill = "Province") +
#   scale_fill_manual(values = province_colors)+
#   theme_minimal()
# 
# a
# 
# ggsave("naive_coi_provinces_ditros.png", a, width = 14, height = 6, bg = "white")
# 
# 
# #naive coi pairwise kruskal wallis provinces
# pairwise_province_naive_coi <- pairwise.wilcox.test(coi_results$naive_coi, 
#                                                     coi_results$province, p.adjust.method = "bonferroni")
# 
# pairwise_province_naive_coi <- melt(pairwise_province_naive_coi[[3]])
# signif_p.pairwise_province_naive_coi<- pairwise_province_naive_coi[pairwise_province_naive_coi$value <0.05 & !is.na(pairwise_province_naive_coi$value),]
# 
# 
# pairwise_province_combinations <- lapply(1:nrow(signif_p.pairwise_province_naive_coi), function(i) {
#   as.character(c(signif_p.pairwise_province_naive_coi[i, "Var1"], 
#                  signif_p.pairwise_province_naive_coi[i, "Var2"]))
# })
# 
# a1 <- ggplot(coi_results, aes(x = province, y = naive_coi, fill = province)) +
#   geom_violin(width = 1, aes(color = region), alpha = 0.4) +
#   geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 11), 
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   ) +
#   scale_fill_discrete(name = "Region") +
#   labs(x = "", y = "Naive COI") +
#   guides(color = FALSE, fill=FALSE)+
#   scale_fill_manual(values = province_colors)+
#   stat_compare_means(comparisons = pairwise_province_combinations, aes(label = after_stat(p.signif)),
#                      method = "wilcox.test")
# 
# a1
# 
# ggsave(paste0("naive_coi_provinces_violin_", SAMPLING,".png"), a1, width = 8, height = 6, bg = "white")
# 
# 
# 
# b <- ggplot(coi_results_region, aes(x = naive_coi, fill = region)) +
#   geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
#   facet_grid( ~ region, scales = "free") +
#   labs(title = "",
#        x = "Naive COI",
#        y = "Frequency",
#        fill = "Region") +
#   theme_minimal() +
#   guides(fill = FALSE)
# 
# b
# 
# ggsave("naive_coi_regions_ditros.png", b, width = 14, height = 6, bg = "white")
# 
# 
# #naive coi pairwise kruskal wallis regions
# pairwise_region_naive_coi <- pairwise.wilcox.test(coi_results_region$naive_coi, 
#                                                   coi_results_region$region, p.adjust.method = "bonferroni")
# 
# pairwise_region_naive_coi <- melt(pairwise_region_naive_coi[[3]])
# signif_p.pairwise_region_naive_coi<- pairwise_region_naive_coi[pairwise_region_naive_coi$value <0.05 & !is.na(pairwise_region_naive_coi$value),]
# 
# 
# pairwise_region_combinations <- lapply(1:nrow(signif_p.pairwise_region_naive_coi), function(i) {
#   as.character(c(signif_p.pairwise_region_naive_coi[i, "Var1"], 
#                  signif_p.pairwise_region_naive_coi[i, "Var2"]))
# })
# 
# 
# b1 <- ggplot(coi_results_region, aes(x = region, y = naive_coi, fill = region)) +
#   geom_violin(width = 1, aes(color = region), alpha = 0.4) +
#   geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 11), 
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   ) +
#   scale_fill_discrete(name = "Region") +
#   labs(x = "", y = "Naive COI") +
#   guides(color = FALSE, fill = FALSE) +
#   stat_compare_means(comparisons = pairwise_region_combinations, aes(label = after_stat(p.signif)),
#                      method = "wilcox.test")
# 
# b1
# 
# ggsave(paste0("naive_coi_regions_violin_", SAMPLING,".png"), b1, width = 8, height = 6, bg = "white")


# #ECOI
#
# c <- ggplot(coi_results, aes(x = post_effective_coi_med, fill = province)) +
#   geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
#   facet_wrap(~ province , scales = "fixed", nrow = 1) + 
#   labs(title = "",
#        x = "Post Effective COI Median",
#        y = "Frequency",
#        fill = "Province") +
#   theme_minimal() +
#   scale_fill_manual(values = province_colors)+
#   guides(fill = FALSE) 
# 
# c
# 
# ggsave("ecoi_provinces_ditros.png", c, width = 14, height = 6, bg = "white")


#ecoi pairwise kruskal wallis provinces
pairwise_province_ecoi <- pairwise.wilcox.test(coi_results$post_effective_coi_med, 
                                               coi_results$province, p.adjust.method = "bonferroni")

pairwise_province_ecoi <- melt(pairwise_province_ecoi[[3]])
signif_p.pairwise_province_ecoi<- pairwise_province_ecoi[pairwise_province_ecoi$value <0.05 & !is.na(pairwise_province_ecoi$value),]

pairwise_province_combinations <- lapply(1:nrow(signif_p.pairwise_province_ecoi), function(i) {
  as.character(c(signif_p.pairwise_province_ecoi[i, "Var1"], 
                 signif_p.pairwise_province_ecoi[i, "Var2"]))
})

c1 <- ggplot(coi_results, aes(x = province, y = post_effective_coi_med, fill = province)) +
  geom_violin(width = 1, aes(color = region), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_discrete(name = "Region") +
  labs(x = "", y = "eCOI") +
  guides(color = FALSE, fill = FALSE) +
  scale_fill_manual(values = province_colors) #+
  #stat_compare_means(comparisons = pairwise_province_combinations, aes(label = after_stat(p.signif)),
  #                   method = "wilcox.test")

#############################3
# Prepare significance labels
signif_p.pairwise_province_ecoi <- signif_p.pairwise_province_ecoi %>%
  mutate(label = case_when(
    value < 0.001 ~ "***",
    value < 0.01 ~ "**",
    value < 0.05 ~ "*",
    TRUE ~ ""
  ))

signif_p.pairwise_province_ecoi <- signif_p.pairwise_province_ecoi[signif_p.pairwise_province_ecoi$value < 0.05,]

# Add significance annotations only for significant results
for (i in 1:nrow(signif_p.pairwise_province_ecoi)) {
  if (signif_p.pairwise_province_ecoi$label[i] != "") {
    c1 <- c1 + geom_signif(
      comparisons = list(c(as.character(signif_p.pairwise_province_ecoi$Var1)[i], as.character(signif_p.pairwise_province_ecoi$Var2)[i])),
      annotations = signif_p.pairwise_province_ecoi$label[i],
      map_signif_level = TRUE,
      y_position = max(coi_results$post_effective_coi_med) + 1 * (i + 1)
    )
  }
}

#############################3


ggsave(paste0("ecoi_provinces_violin_", SAMPLING,".png"), c1, width = 8, height = 7, bg = "white")


# d <- ggplot(coi_results_region, aes(x = post_effective_coi_med, fill = region)) +
#   geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
#   facet_grid( ~ region, scales = "free") +
#   labs(title = "",
#        x = "Post Effective COI Median",
#        y = "Frequency",
#        fill = "Region") +
#   theme_minimal() +
#   guides(fill = FALSE) 
# 
# d
# 
# ggsave("ecoi_regions_ditros.png", d, width = 10, height = 6, bg = "white")


#ecoi pairwise kruskal wallis regions
pairwise_region_ecoi <- pairwise.wilcox.test(coi_results_region$post_effective_coi_med, 
                                             coi_results_region$region, p.adjust.method = "bonferroni")

pairwise_region_ecoi <- melt(pairwise_region_ecoi[[3]])
signif_p.pairwise_region_ecoi<- pairwise_region_ecoi[pairwise_region_ecoi$value <0.05 & !is.na(pairwise_region_ecoi$value),]

pairwise_region_combinations <- lapply(1:nrow(signif_p.pairwise_region_ecoi), function(i) {
  as.character(c(signif_p.pairwise_region_ecoi[i, "Var1"], 
                 signif_p.pairwise_region_ecoi[i, "Var2"]))
})


d1 <- ggplot(coi_results_region, aes(x = region, y = post_effective_coi_med, fill = region)) +
  geom_violin(width = 1, aes(color = region), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_discrete(name = "Region") +
  labs(x = "", y = "eCOI") +
  guides(color = FALSE, fill = FALSE) #+
  #stat_compare_means(comparisons = pairwise_region_combinations, aes(label = after_stat(p.signif)),
  #                   method = "wilcox.test")

#############################3
# Prepare significance labels
signif_p.pairwise_region_ecoi <- signif_p.pairwise_region_ecoi %>%
  mutate(label = case_when(
    value < 0.001 ~ "***",
    value < 0.01 ~ "**",
    value < 0.05 ~ "*",
    TRUE ~ ""
  ))

signif_p.pairwise_region_ecoi <- signif_p.pairwise_region_ecoi[signif_p.pairwise_region_ecoi$value < 0.05,]

# Add significance annotations only for significant results
for (i in 1:nrow(signif_p.pairwise_region_ecoi)) {
  if (signif_p.pairwise_region_ecoi$label[i] != "") {
    d1 <- d1 + geom_signif(
      comparisons = list(c(as.character(signif_p.pairwise_region_ecoi$Var1)[i], as.character(signif_p.pairwise_region_ecoi$Var2)[i])),
      annotations = signif_p.pairwise_region_ecoi$label[i],
      map_signif_level = TRUE,
      y_position = max(coi_results$post_effective_coi_med) + 1 * (i + 1)
    )
  }
}

d1
#############################3

ggsave(paste0("ecoi_coi_pregions_violin_", SAMPLING,".png"), d1, width = 8, height = 7, bg = "white")


# ## %POLYCLONAL
# 
# ##polyclonal percentage
# e <- ggplot(polyclonal_percentage_region, aes(x = region, y = polyclonal_percentage_region, fill = region)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "", y = "% Polyclonal Infections") +
#   #facet_wrap(~region, scales = "fixed", ncol = 3) +
#   #scale_fill_manual(values = c("2022" = "orange")) + 
#   theme_minimal()+
#   guides(fill = FALSE) 
# 
# ggsave(paste0("perc_polyclonal_regions_", SAMPLING, ".png"), e, width = 6, height = 4, bg = "white")
# 
# 
# polyclonal_percentage_province <- merge(polyclonal_percentage_province, unique(coi_results[c("province", "region")]), by="province")
# 
# f <- ggplot(polyclonal_percentage_province, aes(x = province, y = polyclonal_percentage_province,  fill = province))+
#   geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
#   labs(x = "", y = "%Polyclonal Infections") +
#   #facet_wrap(~province, scales = "fixed", nrow = 1) +
#   #scale_fill_manual(values = c("2022" = "orange")) + 
#   theme_minimal()+
#   scale_fill_manual(values = province_colors)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   guides(fill = FALSE) 
# 
# ggsave(paste0("perc_polyclonal_provinces_", SAMPLING, ".png"), f, width = 8, height = 6, bg = "white")


###
# coi results for Simone
coi_for_db <- coi_results[c("NIDA2", "naive_coi", "post_effective_coi_med")]

write.csv(coi_for_db, paste0("coi_for_db_", SAMPLING, ".csv"), row.names = F)

#plot(coi_for_db$naive_coi, coi_for_db$post_effective_coi_med)

scatter_plot <- ggplot(coi_for_db, aes(x = naive_coi, y = post_effective_coi_med)) +
  geom_point(alpha = 0.75) +
  labs(x = "Naive COI", y = "Post Effective COI Median") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("scatter_plot_",SAMPLING,".png"), scatter_plot, width = 8, height = 6, bg = "white")





library(tidyr)
library(purrr)
library(binom)
library(ggsignif)

### POLYCLONAL infections REGIONS

# Preparing data
coi_results_region <- coi_results_region %>%
  mutate(polyclonal = ifelse(post_effective_coi_med > 1.1, 1, 0))

# Function to perform pairwise Fisher's Exact Test
pairwise_fisher_test <- function(data, region1, region2) {
  data1 <- data %>% filter(region == region1)
  data2 <- data %>% filter(region == region2)
  
  table <- matrix(c(
    sum(data1$polyclonal), nrow(data1) - sum(data1$polyclonal),
    sum(data2$polyclonal), nrow(data2) - sum(data2$polyclonal)
  ), nrow = 2, byrow = TRUE)
  
  fisher.test(table)$p.value
}

# Unique regions
regions <- unique(as.character(coi_results_region$region))

# Generate all pairs of regions
region_pairs <- combn(regions, 2, simplify = FALSE)

# Pairwise comparisons
pairwise_results <- tibble(
  region1 = sapply(region_pairs, `[[`, 1),
  region2 = sapply(region_pairs, `[[`, 2),
  p.value = map2_dbl(sapply(region_pairs, `[[`, 1), sapply(region_pairs, `[[`, 2), 
                     ~pairwise_fisher_test(coi_results_region, .x, .y))
)

# Adjust for multiple comparisons (e.g., Bonferroni correction)
pairwise_results <- pairwise_results %>%
  mutate(p.adjusted = p.adjust(p.value, method = "bonferroni"))

pairwise_results


# Calculate proportion of polyclonal infections and confidence intervals
summary_data <- coi_results_region %>%
  mutate(polyclonal = ifelse(post_effective_coi_med > 1.1, 1, 0)) %>%
  group_by(region) %>%
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
summary_data$region <- factor(summary_data$region, levels = c("North", "Centre", "South"))


e<- ggplot(summary_data, aes(x = region, y = prop_poly, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = "",
       x = "",
       y = "Proportion of Polyclonal Infections",
       fill = "Region") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill = FALSE) 
  

# Prepare significance labels
pairwise_results <- pairwise_results %>%
  mutate(label = case_when(
    p.adjusted < 0.001 ~ "***",
    p.adjusted < 0.01 ~ "**",
    p.adjusted < 0.05 ~ "*",
    TRUE ~ ""
  ))
pairwise_results <- pairwise_results[pairwise_results$p.adjusted < 0.05,]

# Add significance annotations only for significant results
for (i in 1:nrow(pairwise_results)) {
  if (pairwise_results$label[i] != "") {
    e <- e + geom_signif(
      comparisons = list(c(pairwise_results$region1[i], pairwise_results$region2[i])),
      annotations = pairwise_results$label[i],
      map_signif_level = TRUE,
      y_position = max(summary_data$prop_poly) + 0.05 * (i + 1)
    )
  }
}

ggsave(paste0("perc_polyclonal_regions_", SAMPLING, ".png"), e, width = 6, height = 5, bg = "white")




### POLYCLONAL infections PROVINCES

# Preparing data
coi_results <- coi_results %>%
  mutate(polyclonal = ifelse(post_effective_coi_med > 1.1, 1, 0))

# Function to perform pairwise Fisher's Exact Test
pairwise_fisher_test <- function(data, region1, region2) {
  data1 <- data %>% filter(province == region1)
  data2 <- data %>% filter(province == region2)
  
  table <- matrix(c(
    sum(data1$polyclonal), nrow(data1) - sum(data1$polyclonal),
    sum(data2$polyclonal), nrow(data2) - sum(data2$polyclonal)
  ), nrow = 2, byrow = TRUE)
  
  fisher.test(table)$p.value
}

# Unique regions
regions <- unique(as.character(coi_results$province))

# Generate all pairs of regions
region_pairs <- combn(regions, 2, simplify = FALSE)

# Pairwise comparisons
pairwise_results <- tibble(
  region1 = sapply(region_pairs, `[[`, 1),
  region2 = sapply(region_pairs, `[[`, 2),
  p.value = map2_dbl(sapply(region_pairs, `[[`, 1), sapply(region_pairs, `[[`, 2), 
                     ~pairwise_fisher_test(coi_results, .x, .y))
)

# Adjust for multiple comparisons (e.g., Bonferroni correction)
pairwise_results <- pairwise_results %>%
  mutate(p.adjusted = p.adjust(p.value, method = "bonferroni"))

pairwise_results

# Calculate proportion of polyclonal infections and confidence intervals
summary_data <- coi_results %>%
  mutate(polyclonal = ifelse(post_effective_coi_med > 1.1, 1, 0)) %>%
  group_by(province) %>%
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
summary_data$province <- factor(summary_data$province, levels = provinces)


f<- ggplot(summary_data, aes(x = province, y = prop_poly, fill = province)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = "",
       x = "",
       y = "Proportion of Polyclonal Infections",
       fill = "Region") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill = FALSE) +
  scale_fill_manual(values = province_colors)


# Prepare significance labels
pairwise_results <- pairwise_results %>%
  mutate(label = case_when(
    p.adjusted < 0.001 ~ "***",
    p.adjusted < 0.01 ~ "**",
    p.adjusted < 0.05 ~ "*",
    TRUE ~ ""
  ))
pairwise_results <- pairwise_results[pairwise_results$p.adjusted < 0.05,]


# Add significance annotations only for significant results
for (i in 1:nrow(pairwise_results)) {
  if (pairwise_results$label[i] != "") {
    f <- f + geom_signif(
      comparisons = list(c(pairwise_results$region1[i], pairwise_results$region2[i])),
      annotations = pairwise_results$label[i],
      map_signif_level = TRUE,
      y_position = max(summary_data$prop_poly) + 0.05 * (i + 1)
    )
  }
}

ggsave(paste0("perc_polyclonal_provinces_", SAMPLING, ".png"), f, width = 8, height = 6.5, bg = "white")

