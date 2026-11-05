# =========================================
# Development of baseline lake whitefish samples for RUBIAS, Moran and Anderson, 2019 (local v0.3)
# Dataset: from the 465-locus baseline (lwfRubias_v1.csv) to the final baseline
# IMPORTANT: (v5.8 in Mich_LakeWhitefish.Rmd): 295 loci, migrants/admixed individuals removed, reporting groups assigned
# IMPORTANT: Each step keeps the Rmd's version name (lwf_rubias_v5, v5.1 ... v5.8) so the two can be compared
# =========================================


library(tidyverse)
library(adegenet)
source("./lwf_functions.R")

#Path to the STRUCTURE binary (mainparams / extraparams / extraparams2 must be in the working directory)
path2structure <- "./structure/STRUCTURE"

# ====

# Input files
CollectionInfo <- read_csv("./baselineIndData.csv")
lwf_rubias_v1 <- read_csv("./lwfRubias_v1.csv") %>%
  mutate(across(5:ncol(.), as.numeric))
RepoGrpMap <- read.csv("./lwf_PopNames.csv")       # Location -> original reporting group (GMZ)
LocSum <- read_delim("./lwf_locus_Summary.csv")    # manual locus quality scores (allele read ratio plots)
LL2Drop <- read_csv("./LinkedLoci_Decision.csv")$Drop %>% unique()   # one locus of each linked pair


# =========================================

# 01 | Select loci to keep and clean collection reporting groups

# Loci that are fixed (<2 alleles) in the baseline samples
LocNames <- colnames(lwf_rubias_v1)[-(1:4)][c(T,F)]
LocusAls <- sapply(LocNames, function(x) length(na.omit(unique(c(lwf_rubias_v1[[x]], lwf_rubias_v1[[paste0(x, ".1")]])))))
# Identify monomorphic loci
FixedLoci <- LocNames[which(LocusAls < 2)] 

LocSum <- LocSum %>%
  mutate(InGPfile = ifelse(LocusName %in% LocNames, "Yes", "No")) %>%
  mutate(Qual_Score_New = ifelse(LocusName %in% FixedLoci, 5, `Qual_Score(1-5)`))

# Retain loci with quality scores 1-2
Loci2Keepv5 <- LocSum %>% 
  filter(Qual_Score_New %in% c(1,2) & InGPfile == "Yes") %>%
  filter(LocusName != "944332_SNP2_124") %>%  # remove new bad locus from HWE issues
  mutate(LocNum = (lapply(str_split(Locus, "_"), "[[", 1) %>% unlist())) %>%
  group_by(LocNum) %>%
  slice_head(n=1) # NOTE: first SNP from each amplicon was taken


# ====

# Annual collections from the same location are pooled (Peshtigo River 2017-2020 included).
# Oconto River (n = 4) is dropped due to low sample size.
lwf_rubias_v5 <- lwf_rubias_v1 %>%
  filter(indiv %in% CollectionInfo$SampleID) %>%
  mutate(collection = plyr::mapvalues(x = indiv,
                                      from = CollectionInfo$SampleID,
                                      to = CollectionInfo$Location, warn_missing = F)) %>%
  filter(!grepl(pattern="Oconto", collection)) %>%
  mutate(repunit = plyr::mapvalues(from = RepoGrpMap$Location,
                                   to = RepoGrpMap$GMZ,
                                   x = collection))

# Replace 0s with NAs
lwf_rubias_v5[lwf_rubias_v5==0] <- NA

# Keep the loci (should be 329), then drop one of each linked pair (should be 295)
Col2Keepv5 <- c(1:4, grep(pattern=paste(Loci2Keepv5$Locus %>% unlist(), collapse="|"), colnames(lwf_rubias_v5)))
lwf_rubias_v5.1 <- lwf_rubias_v5[, Col2Keepv5]

# Identify allele columns belonging to loci designated for removal
LL2Drop_ColIndex <- grep(pattern=paste0("^", LL2Drop, "_", collapse="|"), colnames(lwf_rubias_v5.1))
lwf_rubias_v5.2 <- lwf_rubias_v5.1[, -LL2Drop_ColIndex]

Loci2Keepv5 <- colnames(lwf_rubias_v5.2) %>%
  {.[-(1:4)]} %>%
  {.[c(T,F)]}
length(Loci2Keepv5) # 295


# =========================================

# 02 | STRUCTURE (multiple analyses with different reporting group structure)

## STRUCTURE with all baseline collections, K = 2 (east vs west)

# GOAL: Find west-lineage fish (Q2) sampled in east collections; they are putative migrants
westColls <- c("Fox River", "Riley Bay", "Larsen Reef and Snake Island", "Peshtigo River", "Menominee River",
               "Whitefish Bay", "North and Moonlight Bays", "Big Bay de Noc", "Manistique")
# Treat all remaining collections as eastern lake side for initial assessment
eastColls <- setdiff(unique(lwf_rubias_v5.2$collection), westColls)

# Export all baseline samples to STRUCTURE format
FB <- writeStructureFile(lwf_rubias_v5.2, Loci2Keepv5, "./lwf_21p_295L_Bv5_2_structure.txt")
runStructure(path2structure, "./lwf_21p_295L_Bv5_2_structure.txt", "./Outfiles/lwf_21p_295L",
             K = 2, nInd = FB$nInd, nLoci = length(Loci2Keepv5))

FB_k2_Qm <- readStructureQ("./Outfiles/lwf_21p_295L", K = 2, PopNumbs = FB$PopNumbs) %>%
  alignK2(anchor_collections = westColls, anchor_cluster = "Q2")      # Q2 = west lineage

# Summarize ancestry coefficients across runs
FB_MeanQ <- FB_k2_Qm %>%
  group_by(Sample_ID, Pop, variable_new, Collection) %>%
  summarize(MeanQ = mean(as.numeric(value)),
            SdQ = sd(as.numeric(value)))

# Identify eastern samples assigned to western ancestry
East_AssgQ2 <- FB_MeanQ %>%
  ungroup() %>%
  group_by(Sample_ID, Pop, Collection) %>%
  slice_max(MeanQ) %>%
  filter(Collection %in% eastColls) %>%
  filter(variable_new == "Q2") %>%
  filter(MeanQ >= 0.7)

# Remove putative migrants
lwf_rubias_v5.3 <- lwf_rubias_v5.2 %>%
  filter(!(indiv %in% East_AssgQ2$Sample_ID))

# ====

## STRUCTURE with all west collections, K = 2 (lineages L1 and L2 within GMZ1)

# Export all western samples to STRUCTURE format
NW <- writeStructureFile(lwf_rubias_v5.3, Loci2Keepv5, "./lwf_9p_295L_v5_3_structure.txt", collections = westColls)
runStructure(path2structure, "./lwf_9p_295L_v5_3_structure.txt", "./Outfiles/lwf_9p_295L_v5_3",
             K = 2, nInd = NW$nInd, nLoci = length(Loci2Keepv5))

NW_k2_Qm <- readStructureQ("./Outfiles/lwf_9p_295L_v5_3", K = 2, PopNumbs = NW$PopNumbs) %>%
  alignK2(anchor_collections = c("Whitefish Bay", "North and Moonlight Bays", "Big Bay de Noc"),
          anchor_cluster = "Q1") # Q1 = L1 (northern Green Bay)

# Summarize ancestry coefficients across runs
NW_MeanQ <- NW_k2_Qm %>%
  group_by(Sample_ID, Pop, variable_new, Collection) %>%
  summarize(MeanQ = mean(as.numeric(value)),
            SdQ = sd(as.numeric(value)))

# Assign all western designated fish to its highest Q lineage
NW_Assg <- NW_MeanQ %>%
  ungroup() %>%
  group_by(Sample_ID, Pop, Collection) %>%
  slice_max(MeanQ) %>%
  mutate(collection_new = paste(Collection, " L", gsub(pattern="Q", replacement="", variable_new), sep="")) %>%
  mutate(repunit_new = paste("GMZ1", "L", gsub(pattern="Q", replacement="", variable_new), sep=""))

# Export Q values from STRUCTURE
write_csv(NW_MeanQ %>% ungroup(), "./output/NW_K2_MeanQ.csv")

# Flag western fish with Q < 0.7
NW_Admixed <- NW_Assg %>%
  mutate(MeanQ = round(as.numeric(MeanQ), 2)) %>%
  filter(MeanQ < 0.7)

# Remove admixed wester fish and replace repunit lables
lwf_rubias_v5.3 <- lwf_rubias_v5.3 %>% 
  filter(!(indiv %in% NW_Admixed$Sample_ID)) %>% #remove admixed indiv. from GMZ1
  mutate(collection_new = ifelse(indiv %in% NW_Assg$Sample_ID, plyr::mapvalues(from = NW_Assg$Sample_ID,
                                                                               to = NW_Assg$collection_new,
                                                                               x = indiv), collection)) %>%
  mutate(repunit_new = ifelse(indiv %in% NW_Assg$Sample_ID, plyr::mapvalues(from = NW_Assg$Sample_ID,
                                                                            to = NW_Assg$repunit_new,
                                                                            x = indiv), repunit)) %>%
  relocate(collection_new, .after = collection) %>%
  relocate(repunit_new, .after = repunit) %>%
  select(-c(repunit, collection)) %>%
  rename('repunit' = 'repunit_new',
         'collection' = 'collection_new')

# ====

## STRUCTURE with all east collections, K = 5

E5_run <- 3             # run used for assignment
E5_gmz6_cluster <- "Q2" # cluster that corresponds to the west-shore river collections (GMZ6)

# Export all eastern samples to STRUCTURE format
E <- writeStructureFile(lwf_rubias_v5.3, Loci2Keepv5, "./lwf_12p_295L_v5_3_structure.txt", collections = eastColls)
runStructure(path2structure, "./lwf_12p_295L_v5_3_structure.txt", "./Outfiles/lwf_12p_295L_v5_3",
             K = 5, nInd = E$nInd, nLoci = length(Loci2Keepv5))
E_k5_Qm <- readStructureQ("./Outfiles/lwf_12p_295L_v5_3", K = 5, PopNumbs = E$PopNumbs)

E_Q5_Assgn <- E_k5_Qm %>%
  filter(Run == E5_run) %>%
  group_by(Sample_ID, Pop, Collection) %>%
  slice_max(value)

# Individuals from the GMZ6 collections that are not in the GMZ6 cluster
PutMigIntoGMZ6 <- E_Q5_Assgn %>% 
  filter(variable != E5_gmz6_cluster & (Collection %in% c("Pentwater River", "White River", "Muskegon", "Grand River", "Black River")))
# Individuals not from GMZ6 that are in the GMZ6 cluster
PutMigIntoGMZ2_5 <- E_Q5_Assgn %>% 
  filter(variable == E5_gmz6_cluster & (Collection %in% c("Naubinway", "Epoufette", "Cross Village", "Little Traverse Bay", "Elk Rapids", "Ingalls Point and Ingalls Bay", "Good Harbor Bay")))

# Now remove all those individuals
lwf_rubias_v5.4 <- lwf_rubias_v5.3 %>%
  filter(!(indiv %in% PutMigIntoGMZ2_5$Sample_ID)) %>%
  filter(!(indiv %in% PutMigIntoGMZ6$Sample_ID))

# ====

## STRUCTURE with the seven northeast collections, K = 4 (non-admixture model)

#Define the northeast collections
neColls <- c("Naubinway", "Epoufette", "Cross Village", "Little Traverse Bay", "Elk Rapids", "Ingalls Point and Ingalls Bay", "Good Harbor Bay")
E4_run <- 2
E4_ingalls_cluster <- "Q4"

# Export all samples to STRUCTURE format
NE4 <- writeStructureFile(lwf_rubias_v5.4, Loci2Keepv5, "./lwf_7p_295L_Bv5_4_structure.txt", collections = neColls)
runStructure(path2structure, "./lwf_7p_295L_Bv5_4_structure.txt", "./Outfiles/lwf_7p_295L_Bv5_4",
             K = 4, nInd = NE4$nInd, nLoci = length(Loci2Keepv5), extra = "-e ./extraparams2")
E_k4_Qm <- readStructureQ("./Outfiles/lwf_7p_295L_Bv5_4", K = 4, PopNumbs = NE4$PopNumbs)

# Assign each fish to its highest-Q cluster
E_Q4_Assgn <- E_k4_Qm %>%
  filter(Run == E4_run) %>%
  group_by(Sample_ID, Pop, Collection) %>%
  slice_max(value)

# Flag Ingalls Point/Bay fish
PutMigIntoGMZ5 <- E_Q4_Assgn %>% 
  filter(variable != E4_ingalls_cluster & (Collection %in% c("Ingalls Point and Ingalls Bay")))
PutMigOutGMZ5 <- E_Q4_Assgn %>% 
  filter(variable == E4_ingalls_cluster & (!(Collection %in% c("Ingalls Point and Ingalls Bay"))))

# Remove flagged fish
lwf_rubias_v5.5 <- lwf_rubias_v5.4 %>%
  filter(!(indiv %in% PutMigIntoGMZ5$Sample_ID))

# ====

## STRUCTURE without Ingalls Point/Bay, Elk Rapids, and Good Harbor Bay, K = 3

NE3_run <- 1
NE3_gmz4_cluster <- "Q1"

# Export files and run STRUCTURE K =3
NE3 <- writeStructureFile(lwf_rubias_v5.5, Loci2Keepv5, "./lwf_6p_295L_Bv5_5_structure.txt",
                          collections = setdiff(neColls, "Ingalls Point and Ingalls Bay"))
runStructure(path2structure, "./lwf_6p_295L_Bv5_5_structure.txt", "./Outfiles/lwf_6p_295L_Bv5_5",
             K = 3, nInd = NE3$nInd, nLoci = length(Loci2Keepv5), extra = "-e ./extraparams2")
NE_k3_Qm <- readStructureQ("./Outfiles/lwf_6p_295L_Bv5_5", K = 3, PopNumbs = NE3$PopNumbs)

# Assign each fish to its highest-Q cluster
NE_Q3_Assgn <- NE_k3_Qm %>%
  filter(Run == NE3_run) %>%
  group_by(Sample_ID, Pop, Collection) %>%
  slice_max(value)

# Indinividuals from GMZ4 that aren't in the GMZ4 cluster
PutMigIntoGMZ4 <- NE_Q3_Assgn %>% 
  filter(variable != NE3_gmz4_cluster & (Collection %in% c("Elk Rapids", "Good Harbor Bay")))

lwf_rubias_v5.6 <- lwf_rubias_v5.5 %>%
  filter(!(indiv %in% PutMigIntoGMZ4$Sample_ID))

# Small lineage-specific GMZ1 collections that were not carried forward, remove putative migrants
lwf_rubias_v5.6 <- lwf_rubias_v5.6 %>%
  filter(!(collection %in% c("Larsen Reef and Snake Island L1", "Riley Bay L2", "Riley Bay L1", "Peshtigo River L1")))

# ====

## STRUCTURE with the four remaining northeast collections, K = 2 

# Check only: these collections could not be resolved any further and form one reporting group
neColls2 <- c("Naubinway", "Epoufette", "Cross Village", "Little Traverse Bay")
NE2 <- writeStructureFile(lwf_rubias_v5.6, Loci2Keepv5, "./lwf_4p_295L_Bv5_6_structure.txt", collections = neColls2)
runStructure(path2structure, "./lwf_4p_295L_Bv5_6_structure.txt", "./Outfiles/lwf_4p_295L_Bv5_6",
             K = 2, nInd = NE2$nInd, nLoci = length(Loci2Keepv5), extra = "-e ./extraparams2")


# =========================================

# 03 | Principal Component Analysis (Remove outliers from the northeast GMZ2 collections)

# Genepop file of the current baseline so that adegenet can read it
writeGenepop(lwf_rubias_v5.6, "lwf_22p_295loci_v5.6", "./lwf_22p_295loci_v5.6_GP.gen")
lwf_22Clcts_295loci_genid <- adegenet::read.genepop("./lwf_22p_295loci_v5.6_GP.gen", ncode=2)

PCA_Res <- list()
for (p in neColls2) {
  Pop_inds <- lwf_rubias_v5.6 %>% filter(collection == p) %>% pull(indiv)
  if (length(Pop_inds) < 5) next   # nothing left to evaluate
  Pop_genid <- lwf_22Clcts_295loci_genid[indNames(lwf_22Clcts_295loci_genid) %in% Pop_inds]
  PCA_Res[[p]] <- pcaOutliers(Pop_genid) %>% mutate(Pop = p)
}

# Combine all outlier results
PCA_Outliers <- PCA_Res %>%
  do.call(rbind, .) %>%
  as_tibble()

# Remove outliers
lwf_rubias_v5.7 <- lwf_rubias_v5.6 %>% 
  filter(!indiv %in% (PCA_Outliers %>% filter(Fp < 0.1) %>% select(Ind) %>% unlist())) #remove outliers


# =========================================

# 04 | Export final baseline dataset

lwf_rubias_v5.8 <- lwf_rubias_v5.7 %>%
  filter(!(indiv %in% PutMigOutGMZ5$Sample_ID)) %>%
  mutate(repunit = recode(repunit,
                          'GMZ2' = 'GMZ2_3',
                          'GMZ3' = 'GMZ2_3'))

# Reporting group names as used in the manuscript (L1 = northern lineage, L2 = southern lineage)
lwfRubiasBaseline <- lwf_rubias_v5.8 %>%
  mutate(repunit = recode(repunit,
                          'GMZ1L1' = 'GMZ1N',
                          'GMZ1L2' = 'GMZ1S',
                          'GMZ2_3' = 'GMZ2',
                          'GMZ4' = 'GMZ3',
                          'GMZ5' = 'GMZ4',
                          'GMZ6' = 'GMZ5'))

lwfRubiasBaseline %>% count(repunit)    
write_csv(lwfRubiasBaseline, "./lwfRubiasBaseline.csv")
