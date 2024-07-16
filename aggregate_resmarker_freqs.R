
# currently, allele frequencies for resmarkers are presented as haplotypes. 
# edu needs the frequencies of each resmarker separately, by region and province.
# i need to separate frequencies for each resmarker from the curren haplotypes and sum them up for each pop. easy

library(dplyr)
library(tidyr)

        
#import data
af_2021 <- read.csv("ALLELE_fREQs_DIV_RES_2021.csv")
af_2022_inclures_dry <- read.csv("ALLELE_fREQs_DIV_RES_2022_includes_DRY.csv")

af_all <- rbind(af_2021, af_2022_inclures_dry)


# keep only resmarkers data
resmarkers <- af_all[!grepl("Pf3D7_", af_all$allele), ]
unique(resmarkers$population)

# Split the 'allele' column by "__" into multiple new columns
resmarkers <- resmarkers %>%
  separate(allele, into = c("gene", "position", "variant"), sep = "__", remove = FALSE) %>%
  select(-locus, -allele)

# keep dhps, dhfr
resmarkers <- resmarkers[resmarkers$gene %in% c("dhps", "dhfr"),]

#separate by each variant (by "/")
resmarkers <- resmarkers %>%
  separate_rows(position, variant, sep = "/")

# aggregate each variant by population (sum everyghing)
resmarkers_by_variant <- resmarkers %>%
  group_by(gene, position, variant, population, geo) %>%
  summarize(post_allele_freqs_lower = sum(post_allele_freqs_lower), 
            post_allele_freqs_med = sum(post_allele_freqs_med),
            post_allele_freqs_upper = sum(post_allele_freqs_upper),
            post_allele_freqs_mean = sum(post_allele_freqs_mean))


resmarkers_by_variant <- resmarkers_by_variant %>%
  mutate_if(is.numeric, round, 4)


# dhps 581 comes up as 2 in freq. probs 2 amplicons ugh... just substract 1 and done. they are all monoallelic anyways...

resmarkers_by_variant$post_allele_freqs_mean <- ifelse(resmarkers_by_variant$post_allele_freqs_mean > 1, 1, resmarkers_by_variant$post_allele_freqs_mean)
resmarkers_by_variant$post_allele_freqs_upper <- ifelse(resmarkers_by_variant$post_allele_freqs_upper > 1, 1, resmarkers_by_variant$post_allele_freqs_upper)
resmarkers_by_variant$post_allele_freqs_med <- ifelse(resmarkers_by_variant$post_allele_freqs_med > 1, 1, resmarkers_by_variant$post_allele_freqs_med)
resmarkers_by_variant$post_allele_freqs_lower <- ifelse(resmarkers_by_variant$post_allele_freqs_lower > 1, 1, resmarkers_by_variant$post_allele_freqs_lower)


write.csv(resmarkers_by_variant, "ALLELE_FREQS_RES-F_2021_2022_incluresdry.csv", row.names= F)
