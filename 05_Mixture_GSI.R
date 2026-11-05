# =========================================
# Lake whitefish: RUBIAS genetic stock identification (GSI)
# Toolkit: RUBIAS, Moran and Anderson, 2019 (local v0.3)
# =========================================

library(tidyverse)
library(rubias)
library(coda)
source("./lwf_functions.R")

# =========================================

# 01 | Load baseline data and run the GSI

# Read in files to run rubias (made by 03_Baseline_Development.R and 02_Build_Mixture_Rubias.R)
refDat <- read_csv("./lwfRubiasBaseline.csv")

# Drop the two loci removed from the mixed stock datset in filtering above (line 74 in 02_Build_Mixture_Rubias.R)
refDat2 <- refDat %>%
  select(-c("891923_SNP1_102.1", "891923_SNP1_102", "569374_SNP1_156.1", "569374_SNP1_156")) %>%
  mutate_all(as.character)

# Create flat prior
ncoll <- n_distinct(refDat2$collection)    #Update 'reference'

flat_prior <- tibble(
  collection = "DEFAULT_PI",
  pi_param   = 1 / ncoll)

# ====

# Run GSI
# Mixture collections were analysed grouped by zone + year + season ("SeasonYearZone") and by zone + year ("YearZone")
runMSA <- function(grouping) {
  mixDat <- read_csv(paste0("./output/lwfRubiasMixed_", grouping, ".csv"))
  mixDat2 <- mixDat %>%
    mutate_all(as.character)
  
  ##Run MSA##   six chains
  MCMCReps <- 10000  
  BURNin <- MCMCReps/2
  nChains <- 6
  
  chain_est <- list()
  chain_trace <- list()
  for (chain in 1:nChains) {
    set.seed(chain)
    mix_est <- infer_mixture(reference = refDat2, 
                             mixture = mixDat2,
                             gen_start_col = 5,
                             reps = MCMCReps,
                             burn_in = BURNin,
                             method = "PB",
                             pb_iter=100,
                             pi_prior = flat_prior)  #method PB is for bootstrapping 
    
    ## Analyze MSA ##
    msa <- msaSummary(mix_est, BURNin)
    chain_est[[chain]] <- msa$RG_mix_ests %>% mutate(chain = chain)
    chain_trace[[chain]] <- msa$trace_subset %>% ungroup() %>% mutate(chain = chain)
    
    if (chain == 1) {
      #individual posteriors (from the first chain)
      Test <- mix_est$indiv_posteriors %>%
        select(c(mixture_collection, indiv, repunit, PofZ))
      write_csv(Test, paste0("./output/Individual_posteriors_", grouping, ".csv"))
    }
  }
 


# =========================================

# 02 | Assess chain convergence (Calculate Geman-Rubin potential scale reduction factor for each mixture collection)

# Gelman-Rubin statistic (coda) for each reporting group in each mixture collection
gelman <- bind_rows(chain_trace) %>%
  select(mixture_collection, repunit, sweep, chain, repprop) %>%
  group_by(mixture_collection, repunit) %>%
  group_modify(~{
    chains <- .x %>%
      pivot_wider(names_from = chain, values_from = repprop) %>%
      select(-sweep)
    if (all(apply(chains, 2, sd) == 0)) return(tibble(PSRF = 1))      # reporting group that is never assigned
    tibble(PSRF = gelman.diag(mcmc.list(lapply(chains, as.mcmc)), autoburnin = FALSE)$psrf[1, 1])
  })

gelman %>% ungroup() %>% summarise(max(PSRF))     # should be < 1.1
write_csv(gelman, paste0("./output/GelmanRubin_", grouping, ".csv"))


#############################################################################
#### Pool chains: bootstrap-corrected estimate (mean over chains), 95% CI from the pooled trace ####
RG_mix_ests <- bind_rows(chain_est) %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = mean(repprop),
            bs_corrected_repunit_ppn = mean(bs_corrected_repunit_ppn))

trace_subset <- bind_rows(chain_trace)

BSMod <- msaCI(RG_mix_ests, trace_subset)
write_csv(BSMod, paste0("./output/GSI_estimates_BSMod_", grouping, ".csv"))

##############################################################
#Retrieve Mixture Collection Proportions
SeasonalCollection <- BSMod %>%
  select(mixture_collection, repunit, repprop) %>% 
  mutate(repprop = if_else(repprop < 0,
                           0,
                           repprop)) %>% 
  mutate(repprop = round(repprop, 2)) %>% 
  pivot_wider(names_from = repunit,
              values_from = repprop) 

MC <- SeasonalCollection %>% relocate(mixture_collection,
                                      GMZ1N,
                                      GMZ1S,
                                      GMZ2,
                                      GMZ3,
                                      GMZ4,
                                      GMZ5)
write.csv(MC, paste0("./output/Mixture_Proportions_Collection_", grouping, ".csv"))

####
options(scipen = 999)
Table <- BSMod %>%
  select(mixture_collection, repunit, bs_corrected_repunit_ppn, loCI_BS, hiCI_BS, sd_BS) %>%
  mutate(across(bs_corrected_repunit_ppn:sd_BS, ~round(., 2))) %>%
  filter(bs_corrected_repunit_ppn != 0)
write_csv(Table, paste0("./output/TableS1_GSI_estimates_", grouping, ".csv"))

#################################################################################
######################################P=0#################################
P0 <- trace_subset %>%
  group_by(mixture_collection, repunit) %>%
  mutate(LT_1in1Mill = ifelse(BS_repprop < 0.0000005, 1, 0)) %>%
  summarise(P0 = sum(LT_1in1Mill)/((MCMCReps - BURNin) * nChains))
write_csv(P0, paste0("./output/P0_", grouping, ".csv"))

}


runMSA("SeasonYearZone")
runMSA("YearZone")


