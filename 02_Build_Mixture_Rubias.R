# =========================================
# Data preparation and create of mixed-stock genotypes for RUBIAS, Moran and Andersone, 2019 (local v0.3)
# Dataset: Lake whitefish (Coregonus clupeaformis) genepop files from 01_GTscore_genotyping.R
# =========================================

library(readxl)
library(poppr)
library(adegenet)
library(tidyverse)
library(radiator)
# =========================================

# 01 | Import and harmonize repeated sequencing runs

# Read each genepop file
dat9.genind <- read.genepop("./genepop/2204_Nov_REDO_Cclu_R1_polyGenResults_singleSNP_genepop.gen",
                            ncode = 2,
                            quiet = F) 

dat10.genind <- read.genepop("./genepop/2204_Sep_REDO_Cclu_OLD_R1_polyGenResults_singleSNP_genepop_1.gen",
                             ncode = 2,
                             quiet = F)

# Identify loci present in both repeated sequencing runs
common_loci <- intersect(locNames(dat9.genind), locNames(dat10.genind))

# Subset both genind objects to only include common loci
dat9_filtered.genind <- dat9.genind[loc = common_loci]
dat10_filtered.genind <- dat10.genind[loc = common_loci]

# Retain shared loci
datREDO.genind <- repool(dat9_filtered.genind, dat10_filtered.genind)
# =========================================


# 02 | Harmonize repeat-run genotypes with primary genepop file

# Read genepop file
dat8.genind <- read.genepop("./genepop/Cclu_All_polyGenResults_singleSNP_genepop.gen",
                            ncode = 2,
                            quiet = F)

# Identify loci present in primary genepop file and pooled repeat-run dataset
common_loci <- intersect(locNames(dat8.genind), locNames(datREDO.genind))

# Subset both genind objects to only include common loci
dat8_filtered.genind <- dat8.genind[loc = common_loci]
datREDO_filtered.genind <- datREDO.genind[loc = common_loci]


# ====

# Pull out sample identifiers in both datasets
samplesDat8 <- indNames(dat8_filtered.genind)
samplesREDO <- indNames(datREDO_filtered.genind)

# Identify duplicate samples
duplicates <- intersect(samplesDat8, samplesREDO)

# Remove overlapping samples from genind8 only
genind_filtered_final <- dat8_filtered.genind[!indNames(dat8_filtered.genind) %in% duplicates,]

# Combine files together to create final retained genind object
dat.genind <- repool(genind_filtered_final,
                     datREDO_filtered.genind)


# =========================================

# 03 | Missing-data quality control and duplicate sample review

# Exclude individuals with >= 30% missing genotyping calls
# Exclude loci with >= 30% missing genotype calls
dat.genind_missingREM <- dat.genind %>%
  missingno(type = "genotype",
            cutoff = 0.3) %>% 
  missingno(type = "loci",
            cutoff = 0.3)

# Inspect samples for duplicates after pooling
removeSamples <- tibble(value = indNames(dat.genind_missingREM)) %>%
  filter(duplicated(value) | duplicated(value, fromLast = TRUE)) %>%
  distinct() %>%
  arrange(value)

if (nrow(removeSamples) > 0) {
  stop("Duplicate individual identifiers remain after pooling. Review removeSamples before continuing.")
}

# Create final genind object post filtering and duplicate sample assessment
dat.genind_duplicatesRemoved <- dat.genind_missingREM


# =========================================

# 04 | Join sample metadata and select mixed-stock samples

# Import sample metadata
metaDat <- read_excel("./sampleInfoKey_edited.xlsx") %>% 
  filter(SampleID %in% indNames(dat.genind_duplicatesRemoved)) 

# Retain metadata for genotyped individuals and retain indivduals in genind that are in metadata
dat.genind_metaDat <- dat.genind_duplicatesRemoved[indNames(dat.genind_duplicatesRemoved) %in% metaDat$SampleID]

# Align metadata with genind file
metaDat <- metaDat %>% 
  arrange(match(SampleID, rownames(dat.genind_metaDat@tab)))

# Designate mixed stock sample type as population factor in genind
dat.genind_metaDat@pop <- as.factor(metaDat$SampleType)
# Retain only mixed stock samples (MSA)
dat.mixed <- dat.genind_metaDat[pop = "MSA"]


# =========================================

# 05 | Restrict mixture dataset to two-allele-column loci (rubias file creation did no work with loci represented by one allele column)

# Count allele columns per locus in genind file
allele_counts <- table(adegenet::locFac(dat.mixed))
# Identify loci with two allele columns
keep_loci <- names(allele_counts[allele_counts == 2])
# Match locus names to their positions in the genind object
keep_idx <- which(locNames(dat.mixed) %in% keep_loci)
# Remove loci with only one allele column
dat.mixed.filtered <- dat.mixed[, loc = keep_idx]


# =========================================

# 06 | Retain metadata for MSA samples and align it with genotype rows for RUBIAS file creation

# Subset metadata identifiers
Dat.mixed <- metaDat %>% 
  filter(SampleType %in% c("MSA")) %>% 
  arrange(match(SampleID, rownames(dat.mixed.filtered@tab))) %>%  
  mutate(MonthCollected = replace_na(MonthCollected, "NA")) %>% 
  mutate(msaSeason = case_when(
    MonthCollected %in% c("Feb", "Mar", "Apr", "May") ~ "Spring",
    MonthCollected %in% c("Jun", "Jul", "Aug") ~ "Summer",
    MonthCollected %in% c("Sep", "Oct") ~ "Fall",
    MonthCollected == "NA" ~ Season)) %>% 
  mutate(Year = if_else(Year %in% c(1993, 1997), "1990s", as.character(Year))) %>%  # two 1990s collections pooled
  mutate(msaLotID = paste(Location,        ### Make changes here to alter groupings
                          Year,            ### (e.g. drop Season to pool by zone and year)
                          msaSeason))


# =========================================

# 07 | Convert genotypes into paired allele columns for RUBIAS, maintaining associated metadata

# Covert genind object into table with individual marker calls
tmp <- adegenet::genind2df(dat.mixed.filtered, sep = "") %>%
  as_tibble(rownames = "INDIVIDUALS")

# Arrange table where each row is an individual sample
Dat.mixed2 <- Dat.mixed %>% 
  arrange(match(SampleID, tmp$INDIVIDUALS)) 

# Specify allele columns
tmp2 <- tmp %>% 
  select(-INDIVIDUALS, -pop)

# Separate genotypes into two-character alleles
outTib <- tibble(a = 1:nrow(tmp2)) %>% 
  select()

for (i in 1:ncol(tmp2)) {
  out <- tmp2 %>% 
    select(i) %>% 
    separate(1,
             into = c(colnames(tmp2)[i], paste0(colnames(tmp2)[i], ".1")),
             sep = 2) 
  
  outTib <- bind_cols(outTib, out)
}

#Add sample metadata to allele columns
lwfMixed <- outTib %>% 
  mutate(across(everything(), ~na_if(., "000"))) %>% 
  mutate(sample_type = "mixture", .before = 1) %>% 
  mutate(repunit = NA, .after = 1) %>% 
  mutate(collection = Dat.mixed2$msaLotID, .after = 2) %>% 
  mutate(indiv = Dat.mixed2$SampleID, .after = 3) %>% 
  mutate(repunit = as.character(repunit))

# Output constructed RUBIAS file
write.csv(lwfMixed, "./MixtureSamples.csv")

