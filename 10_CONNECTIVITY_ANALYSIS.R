
library(haven)
library(dcifer)
library(dplyr)
library(ggplot2)
library(reshape2)

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

#######################################################
# 12.- IBD: Proportion of related pairwise infections using IBD between provinces and regions
#######################################################

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


## CONNECTIVITY ANALYSIS
pardef <- par(no.readonly = TRUE)

## SAMPLING samples ##
#format data
dsmp <- formatDat(combined_df_merged, svar = "NIDA2", lvar = "locus", avar = "allele") #CHANGED pseudo_cigar for allele. DIFF RESULTS?
str(dsmp, list.len = 2)

# format metadata
meta <- unique(combined_df_merged[c("NIDA2", "region", "province")])
meta <- meta[match(names(dsmp), meta$NIDA2), ]  # order samples as in dsmp

#estimate naive coi
lrank <- 2
coi   <- getCOI(dsmp, lrank = lrank)
min(coi)

#estimate allele freqs
afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 
str(afreq, list.len = 2)

#order provinces from north to wouth
nsite     <- table(meta$province)[provinces]
ord       <- order(factor(meta$province, levels = provinces))
dsmp <- dsmp[ord]
coi  <- coi[ ord]

#calculate ibd
decifer_results <- ibdDat(dsmp, coi, afreq,  pval = TRUE, confint = TRUE, rnull = 0, 
                     alpha = 0.05, nr = 1e3)  

saveRDS(decifer_results, paste0("decifer_results_", SAMPLING, ".RDS")) ## AQUÍ VOY
decifer_results <- readRDS(paste0("decifer_results_", SAMPLING, ".RDS"))

pdf(paste0("decifer_results_plot", SAMPLING, ".pdf"), width = 15, height = 15) 

layout(matrix(1:2, 1), width = c(15, 1))
par(mar = c(1, 1, 2, 1))
alpha <- 0.01         
nsmp  <- length(dsmp)
atsep <- cumsum(nsite)[-length(nsite)]
isig  <- which(decifer_results[, , "p_value"] <= alpha, arr.ind = TRUE)
dmat  <- decifer_results[, , "estimate"]
dmat[upper.tri(dmat)] <- t(dmat)[upper.tri(t(dmat))] 

plotRel(dmat, isig = isig, draw_diag = TRUE, alpha = alpha, idlab = FALSE, side_id = c(2, 3), srt_id = c(25, 65), lwd_diag = 0.5, border_sig = "darkviolet")

abline(v = atsep, h = atsep, col = "gray45", lty = 5)
atclin <- cumsum(nsite) - nsite/2
mtext(provinces, side = 3, at = atclin, line = 0.2)
mtext(provinces, side = 2, at = atclin, line = 0.2)

par(mar = c(1, 0, 2, 3))
plotColorbar()

dev.off()


### EXAMINE PAIRS OF SAMPLES!
n_samples <- dim(decifer_results)[1]

#extract info
estimates_df <-as.data.frame(decifer_results[1:n_samples, 1:n_samples, "estimate"])
estimates_p <-as.data.frame(decifer_results[1:n_samples, 1:n_samples, "p_value"])
estimates_CI_lower <- as.data.frame(decifer_results[1:n_samples, 1:n_samples, "CI_lower"])
estimates_CI_upper <- as.data.frame(decifer_results[1:n_samples, 1:n_samples, "CI_upper"])


#sort info into a single df
estimates_df_long <- reshape2::melt(as.matrix(estimates_df))
names(estimates_df_long) <- c("sample1", "sample2", "estimate")
estimates_p_long <- reshape2::melt(as.matrix(estimates_p))
names(estimates_p_long) <- c("sample1", "sample2", "p_value")
estimates_CI_lower_long <- reshape2::melt(as.matrix(estimates_CI_lower))
names(estimates_CI_lower_long) <- c("sample1", "sample2", "CI_lower")
estimates_CI_upper_long <- reshape2::melt(as.matrix(estimates_CI_upper))
names(estimates_CI_upper_long) <- c("sample1", "sample2", "CI_upper")

merged_df <- merge(estimates_df_long, estimates_p_long, by = c("sample1", "sample2"))
merged_df <- merge(merged_df, estimates_CI_lower_long, by = c("sample1", "sample2"))
merged_df <- merge(merged_df, estimates_CI_upper_long, by = c("sample1", "sample2"))


saveRDS(merged_df, paste0("decifer_results_TABLE_", SAMPLING, ".RDS"))
merged_df <- readRDS(paste0("decifer_results_TABLE_", SAMPLING, ".RDS"))

#remove NA rows
merged_df <- merged_df[complete.cases(merged_df),]

#keep significant cases
merged_df_signif <- merged_df[merged_df$p_value < 0.01,]
merged_df_signif <- merged_df[merged_df$estimate > 0.25,]

#merge with geo locations for each pair of samples compared
merged_df_signif_geo <- merge(merged_df_signif, combined_df_merged[, c("NIDA2", "province", "region", "VOC_436_581")], by.x = "sample1", by.y = "NIDA2")
colnames(merged_df_signif_geo)[7:9] <- c("province_s1", "region_s1", "VOC_s1")
merged_df_signif_geo<- distinct(merged_df_signif_geo)

merged_df_signif_geo <- merge(merged_df_signif_geo, combined_df_merged[, c("NIDA2", "province", "region", "VOC_436_581")], by.x = "sample2", by.y = "NIDA2")
colnames(merged_df_signif_geo)[10:12] <- c("province_s2", "region_s2", "VOC_s2")
merged_df_signif_geo<- distinct(merged_df_signif_geo)

combos_iniciales<- paste0(merged_df_signif[,1], merged_df_signif[,2])
combos_finales<- paste0(merged_df_signif_geo[,1], merged_df_signif_geo[,2])

#sanity check
if (length(combos_finales %in% combos_finales) == dim(merged_df_signif_geo)[1]){
  print("merge was successful.")
} else{
  print("grab a coffee.")
}

#make connectivity columns
merged_df_signif_geo$conn_provinces <- paste0(merged_df_signif_geo$province_s1, "_", merged_df_signif_geo$province_s2)
merged_df_signif_geo$conn_regions <- paste0(merged_df_signif_geo$region_s1, "_", merged_df_signif_geo$region_s2)

# Sort the data frame by estimate in descending order
sorted_df <- merged_df_signif_geo %>% arrange(desc(estimate))


#remove Dry for connectivity comparisons
sorted_df<- sorted_df[!grepl("Dry", sorted_df$conn_provinces), ] #REMOVE TO INCLUDE DRY SEASON IN THE CONNECTIVITY COMPARISON

#write.csv( sorted_df, "connectivity_2022.csv", row.names = F)

## PROVINCES CONNECTIVITY

# Calculate the median for each group
median_data <- aggregate(estimate ~ conn_provinces, sorted_df, median)

# Reorder conn_regions based on the median values in descending order
sorted_df$conn_provinces <- factor(sorted_df$conn_provinces, levels = median_data[order(-median_data$estimate), "conn_provinces"])

# Plot with legend and sorted x-axis
prov_conn<- ggplot(sorted_df, aes(x = conn_provinces, y = estimate, fill = conn_regions)) +
  geom_violin(width = 1, aes(color = conn_regions), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = conn_regions), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  scale_fill_discrete(name = "Province") +  # Customize legend title
  ggtitle("Province Connectivity") +
  xlab("")+
  ylab("IBD")+
  guides(color =FALSE)

prov_conn

ggsave(paste0("province_connectivity_", SAMPLING, ".png"), prov_conn, width = 14, height = 7, bg = "white")

## REGIONS CONNECTIVITY

# Calculate the median for each group
median_data <- aggregate(estimate ~ conn_regions, sorted_df, median)

# Reorder conn_regions based on the median values in descending order
sorted_df$conn_regions <- factor(sorted_df$conn_regions, levels = median_data[order(-median_data$estimate), "conn_regions"])

# Plot with legend and sorted x-axis
reg_conn <- ggplot(sorted_df, aes(x = conn_regions, y = estimate, fill = conn_regions)) +
  geom_violin(width = 1, aes(color = conn_regions), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = conn_regions), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("Region Connectivity") +
  xlab("")+
  ylab("IBD")+
  guides(color = FALSE, fill =FALSE) 

reg_conn

ggsave(paste0("region_connectivity_", SAMPLING, ".png"), reg_conn, width = 12, height = 8, bg = "white")

# pairwise proportions of related infections #

#number of pairwise significantly related (IBD) infections for provinces and regions
table(sorted_df$conn_provinces)
table(sorted_df$conn_regions)


#sample sizes
sample_size_provinces <- combined_df_merged %>%
  group_by(province) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_provinces

combined_df_merged_nodry <- combined_df_merged %>%
  filter(!grepl("Dry", province))

sample_size_regions <- combined_df_merged_nodry %>%
  group_by(year, region) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_regions


#calculate percentages for provinces
ibd_samples_provinces <- sorted_df %>%
  group_by(paste0(province_s1, "_", province_s2)) %>%
  summarize(ibd_samples = length(unique(c(sample1, sample2))),
            province_s1 = province_s1,
            province_s2 = province_s2) %>%
  distinct()

ibd_samples_provinces <- ibd_samples_provinces[,-1]

ibd_samples_provinces <- merge(ibd_samples_provinces, sample_size_provinces, by.x = "province_s1", by.y = "province")
colnames(ibd_samples_provinces)[4] <- "province_s1_ss"

ibd_samples_provinces <- merge(ibd_samples_provinces, sample_size_provinces, by.x = "province_s2", by.y = "province")
colnames(ibd_samples_provinces)[5] <- "province_s2_ss"

ibd_samples_provinces$total_samples_pairwise <- ifelse(ibd_samples_provinces$province_s1 != ibd_samples_provinces$province_s2, rowSums(ibd_samples_provinces[c("province_s1_ss", "province_s2_ss")]), ibd_samples_provinces$province_s1_ss)

ibd_samples_provinces$perc_ibd_samples_pairwise <- ibd_samples_provinces$ibd_samples / ibd_samples_provinces$total_samples_pairwise

ibd_samples_provinces$pairwise_comparison <- paste0(ibd_samples_provinces$province_s1, "_", ibd_samples_provinces$province_s2)

# Assuming ibd_samples_regions is your data frame
ibd_samples_provinces <- ibd_samples_provinces %>%
  arrange(desc(perc_ibd_samples_pairwise))  # Sort the data frame by perc_ibd_samples_pairwise in descending order

# Reorder the levels of the pairwise_comparison factor
ibd_samples_provinces$pairwise_comparison <- factor(ibd_samples_provinces$pairwise_comparison, 
                                                    levels = ibd_samples_provinces$pairwise_comparison)

ibd_samples_provinces <- merge(ibd_samples_provinces, sorted_df[c("conn_provinces", "conn_regions")], by.x = "pairwise_comparison", by.y = "conn_provinces")

ibd_samples_provinces <- distinct(ibd_samples_provinces)

# Create the bar plot
prop_ibd_prov <- ggplot(ibd_samples_provinces, aes(x = pairwise_comparison, y = perc_ibd_samples_pairwise, fill = conn_regions)) +
  geom_bar(stat = "identity", alpha =0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "proportion of samples") +
  #ggtitle("proportion of significantly related IBD samples")+
  guides(color = FALSE) 

prop_ibd_prov

ggsave(paste0("prov_prop_IBD_samples_", SAMPLING, ".png"), prop_ibd_prov, width = 16, height = 10, bg = "white")

#calculate percentages for regions
ibd_samples_regions <- sorted_df %>%
  group_by(paste0(region_s1, "_", region_s2)) %>%
  summarize(ibd_samples = length(unique(c(sample1, sample2))),
            region_s1 = region_s1,
            region_s2 = region_s2) %>%
  distinct()

ibd_samples_regions <- ibd_samples_regions[,-1]

ibd_samples_regions <- merge(ibd_samples_regions, sample_size_regions, by.x = "region_s1", by.y = "region")
ibd_samples_regions <- ibd_samples_regions[,-4]
colnames(ibd_samples_regions)[4] <- "region_s1_ss"

ibd_samples_regions <- merge(ibd_samples_regions, sample_size_regions, by.x = "region_s2", by.y = "region")
ibd_samples_regions <- ibd_samples_regions[,-5]
colnames(ibd_samples_regions)[5] <- "region_s2_ss"

ibd_samples_regions$total_samples_pairwise <- ifelse(ibd_samples_regions$region_s1 != ibd_samples_regions$region_s2, rowSums(ibd_samples_regions[c("region_s1_ss", "region_s2_ss")]), ibd_samples_regions$region_s1_ss)

ibd_samples_regions$perc_ibd_samples_pairwise <- ibd_samples_regions$ibd_samples / ibd_samples_regions$total_samples_pairwise

ibd_samples_regions$pairwise_comparison <- paste0(ibd_samples_regions$region_s1, "_", ibd_samples_regions$region_s2)

# Assuming ibd_samples_regions is your data frame
ibd_samples_regions <- ibd_samples_regions %>%
  arrange(desc(perc_ibd_samples_pairwise))  # Sort the data frame by perc_ibd_samples_pairwise in descending order

# Reorder the levels of the pairwise_comparison factor
ibd_samples_regions$pairwise_comparison <- factor(ibd_samples_regions$pairwise_comparison, 
                                                  levels = ibd_samples_regions$pairwise_comparison)

# Create the bar plot
prop_ibd_reg <- ggplot(ibd_samples_regions, aes(x = pairwise_comparison, y = perc_ibd_samples_pairwise, fill = pairwise_comparison)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Pairwise Comparison", y = "Percentage of IBD Samples") +
  ggtitle("Percentage of IBD Samples by Pairwise Comparison")+
  guides(color = FALSE, fill =FALSE) 

prop_ibd_reg

ggsave(paste0("region_prop_IBD_samples_", SAMPLING, ".png") , prop_ibd_reg, width = 16, height = 10, bg = "white")


#pairwise IBD between samples with variants of concern and wildtype parasites from the same or different areas.
#remove "no_genotype" category
sorted_df_full_geno <- sorted_df[sorted_df$VOC_s1 != "no_genotype",]
sorted_df_full_geno <- sorted_df_full_geno[sorted_df_full_geno$VOC_s2 != "no_genotype",]

#mean IBD for each pair of genotypes for each pairwise comparison
sorted_df_full_geno_summary <- sorted_df_full_geno %>% 
  group_by(conn_provinces, pairwise_geno = paste0(VOC_s1, "_", VOC_s2)) %>% 
  summarize(IBD = estimate) %>%
  mutate(geno = ifelse(grepl("WT_WT", pairwise_geno), "WT_vs_WT", "WT_vs_res"))

# 
# # Calculate the median IBD for each conn_provinces value
# median_IBD <- aggregate(IBD ~ conn_provinces, data = sorted_df_full_geno_summary[sorted_df_full_geno_summary$geno == "WT_vs_res",], FUN = median)
# 
# # Reorder the levels of conn_provinces based on the median IBD of WT_vs_WT category
# sorted_df_full_geno_summary$conn_provinces <- factor(sorted_df_full_geno_summary$conn_provinces, 
#                                                      levels = median_IBD[order(-median_IBD$IBD), "conn_provinces"])

# Plot with conn_provinces sorted by median IBD of WT_vs_WT category
res_wt <- ggplot(sorted_df_full_geno_summary, aes(x = conn_provinces, y = IBD, fill = geno)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "IBD") +
  ggtitle("")

res_wt

ggsave(paste0("province_IBD_res_vs_WT_", SAMPLING, ".png"), res_wt, width = 24, height = 12, bg = "white")


# sorted_df_full_geno_summary_WT_only<- sorted_df_full_geno_summary %>%
#   filter(grepl("Sofala", conn_provinces))

ggplot(sorted_df_full_geno_summary, aes(x = geno, y = IBD, fill = pairwise_geno)) +
  geom_boxplot() +
  facet_wrap(~ conn_provinces, scales = "free", nrow = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "IBD") +
  ggtitle("IBD by Connection Provinces and Pairwise Genotype")+
  guides(color = FALSE) 

