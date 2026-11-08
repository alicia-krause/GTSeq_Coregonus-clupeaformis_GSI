# =========================================
# Lake whitefish sample size analysis
# GOAL: Sample 2022 mixed-stock collections from three zones and simulate different mixture sizes to assess sensitivity of RUBIAS GSI analyses
# =========================================

library(tidyverse)
library(rubias)
source("./lwf_functions.R")

# =========================================

# 01 | Load baseline and mixed stock datasets

refDat <- read_csv("./lwfRubiasBaseline.csv") %>%
  select(-c("891923_SNP1_102.1", "891923_SNP1_102", "569374_SNP1_156.1", "569374_SNP1_156")) 

mixDat <- read_csv("./output/lwfRubiasMixed_SeasonYearZone.csv") %>% 
  select(colnames(refDat))

# Make all columns characters
refDat <- refDat %>% mutate_all(as.character)
mixDat <- mixDat %>% mutate_all(as.character)

# Define GMZ reporting group order and coloration
RepGroups <- c("GMZ1N", "GMZ1S", "GMZ2", "GMZ3", "GMZ4", "GMZ5")
RepColors <- c("#A9D18E", "darkgreen", "#C48AC1", "#9DC3E6", "#DBAF5C", "#ADB9CA")

# Create flat prior
ncoll <- n_distinct(refDat$collection)    #Update 'reference'

flat_prior <- tibble(
  collection = "DEFAULT_PI",
  pi_param   = 1 / ncoll)

# =========================================

# 02 | Empirical analysis mixed-stock collection dataset generation

# Sample sets: All 2022 samples (all seasons pooled) from WI-01, WFM-01 and WFM-02 (n = 162-296).
SSzones <- c("WFM-01", "WI-01", "WI-02")
# Set sampling incrments, and MCMC parameters
SSincrement <- 25
SSreps <- 100
MCMCReps <- 10000
BURNin <- MCMCReps/2

# Create dataframe of 2022 samples, splitting zone/year/season collection labels
SSdat <- mixDat %>%
  separate(collection, c("MZ", "Year", "Season"), sep = " ", remove = FALSE) %>%
  filter(Year == 2022, MZ %in% SSzones)

# Confirm sample sizes
SSdat %>% count(MZ) 


# =========================================

# 03 | Empirical size sampling and RUBIAS mixture estimations

# GOAL: For each MZ, gather varying sample size groupings (25 increments), 100 replicates for each partial sampling and 1 full run with all samples
  # All replicates are run with RUBIAS and GMZ reporting group assignment proportions/bootstrapped CI are summarized
SSres <- list() 
for (zone in SSzones) {
  zoneDat <- SSdat %>% filter(MZ == zone) %>% mutate(collection = zone) # Combine season with zone labels
  nZone <- nrow(zoneDat)
  
  for (n in unique(c(seq(SSincrement, nZone, by = SSincrement), nZone))) {
    # all samples: nothing to resample
    nrep <- ifelse(n == nZone, 1, SSreps)
    for (rep in 1:nrep) {
      set.seed(n*1000 + rep) # Make subsamples reproducable
      sub <- zoneDat[sample(1:nZone, n, replace = F), ] %>% select(colnames(refDat))
      sub$collection <- zone
      
      mix_est <- infer_mixture(reference = refDat, # Estimate mixture proportions
                               mixture = sub,
                               gen_start_col = 5,
                               reps = MCMCReps,
                               burn_in = BURNin,
                               method = "PB",
                               pb_iter = 100,
                               pi_prior = flat_prior)
      
      msa <- msaSummary(mix_est, BURNin) # Report reporting group estimates and bootstrap uncertainty
      SSres[[length(SSres) + 1]] <- msaCI(msa$RG_mix_ests, msa$trace_subset) %>%
        ungroup() %>%
        mutate(zone = zone, n = n, rep = rep) 
    }
  }
}

SSres <- bind_rows(SSres) # Save subsampling results
write_csv(SSres, "./output/SampleSize_subsamples.csv")


# =========================================

# 04 | Summarize empirical sample size summary estimates and CI

# Mean estimate and mean 95% CI across the replicates; lines mark +/-10% of the mean estimate at the largest n
SSsum <- SSres %>%
  group_by(zone, repunit, n) %>%
  summarise(mean_est = mean(bs_corrected_repunit_ppn),
            mean_loCI = mean(loCI_BS),
            mean_hiCI = mean(hiCI_BS),
            sd_est = sd(bs_corrected_repunit_ppn)) %>%
  group_by(zone, repunit) %>%
  mutate(ref_mean = mean_est[n == max(n)],
         band_lo = ref_mean * 0.9,
         band_hi = ref_mean * 1.1)
write_csv(SSsum, "./output/SampleSize_summary.csv")

#Generate plot showing average estimated proportions, CI limits, and +/-10% of the mean estimate at the largest n
SSsum %>%
  ggplot(aes(x = n, y = mean_est, color = repunit)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_loCI, ymax = mean_hiCI, fill = repunit), alpha = 0.2, color = NA) +
  geom_hline(aes(yintercept = band_lo, color = repunit), linetype = "dotted") +
  geom_hline(aes(yintercept = band_hi, color = repunit), linetype = "dotted") +
  scale_color_manual(values = RepColors) +
  scale_fill_manual(values = RepColors) +
  facet_grid(repunit ~ zone, scales = "free_y") +
  xlab("Sample size") +
  ylab("Bootstrap-corrected mixture proportion") +
  theme_bw(base_size = 15)

ggsave("./output/Fig6_SampleSize.png", width = 9, height = 11, dpi = 300) #Save plot


# =========================================

# 05 | Simulated leave-one-out analysis over a range of mixture sizes

# Equal Dirichlet parameter for all six reporting groups (alpha_repunit = 1), alpha_collection = 10
SSsizes <- c(25, 50, 75, 100, 150, 200, 300)

# Simulate 100 mixture runs of the selected sample size
loo <- list()
for (SS in SSsizes) {
  Results_lwf_sims <- assess_reference_loo(reference = refDat, 
                                           gen_start_col = 5, 
                                           reps = 100, 
                                           mixsize = SS,
                                           alpha_repunit = 1,
                                           alpha_collection = 10)
  loo[[as.character(SS)]] <- Results_lwf_sims %>% # Group by reporting group GMZ
    mutate(repunit_f = factor(x = repunit, levels = unlist(RepGroups))) %>%
    group_by(repunit_scenario, iter, repunit_f) %>% 
    summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
    mutate(repu_n_prop = repu_n / sum(repu_n)) %>%
    mutate(mixsize = SS)
}
SS_all <- bind_rows(loo) # Comgine all simulation results


# =========================================

# 06 | Calculate root mean squared error (RMSE) of the estimated proportions

# GOAL: Assess differences between true and inferred mixture proportions
RM_Test <- SS_all %>%
  group_by(mixsize, repunit_f) %>%
  summarise(RMSE = sqrt(mean((true_repprop - repprop_posterior_mean)^2)))
write_csv(RM_Test, "./output/SampleSize_RMSE.csv")

# Plot RMSE by mixture size
RM_Test %>%
  ggplot(aes(x = mixsize, y = RMSE, color = repunit_f)) +
  geom_line(lwd = 1) +
  scale_color_manual(values = RepColors) +
  scale_x_continuous(breaks = seq(0, 300, by=50), expand = c(0,0)) +
  ylab("Root Mean Squared Error (RMSE)") +
  xlab("Sample Size") +
  labs(color = "Reporting Unit") +
  theme_bw(base_size = 15)

ggsave("./output/FigS7_RMSE.png", width = 7, height = 5, dpi = 300) # Save Plot
