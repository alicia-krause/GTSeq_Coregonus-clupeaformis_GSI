# =========================================
# Lake whitefish: RUBIAS baseline evaluation
# Analyses: 1) self-assignment, 2) 100% single-reporting group mixture simulations, and 3) simulated fishery mixtures
# Toolkit: RUBIAS, Moran and Anderson, 2019 (local v0.3)
# =========================================

library(tidyverse)
library(rubias)
library(MCMCpack)
dir.create("./output", showWarnings = FALSE)

# =========================================

# 01 | Configuring genetic management zone (GMZ) plotting

# Read in baseline data
refDat <- read_csv("./lwfRubiasBaseline.csv")
refDat[] <- lapply(refDat, as.character)

# Define GMZ reporting group order and coloration
RepGroups <- c("GMZ1N", "GMZ1S", "GMZ2", "GMZ3", "GMZ4", "GMZ5")
RepColors <- c("#A9D18E", "darkgreen", "#C48AC1", "#9DC3E6", "#DBAF5C", "#ADB9CA")


# =========================================

# 02 | Individual self-assignment (leave-on-out analysis)

# Run RUBIAS self-assignment
sa_lwf <- self_assign(reference = refDat, 
                      gen_start_col = 5)

# Create likelihood table from results
sa_to_repu <- sa_lwf %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarize(repu_scaled_like = sum(scaled_likelihood))

# Identify most likely reporting group GMZ for each sample
sa_to_repu.maxLikelihood <-
  sa_to_repu %>% 
  group_by(indiv) %>% 
  filter(repu_scaled_like == max(repu_scaled_like))

# Make confusion matrix
SA_rates <- sa_to_repu.maxLikelihood %>% 
  ungroup() %>% 
  select(indiv, repunit, inferred_repunit) %>%
  group_by(repunit) %>% 
  mutate(total=n()) %>% 
  group_by(repunit, inferred_repunit) %>% 
  mutate(prop=n()/total) %>% 
  ungroup() %>% 
  select(repunit, inferred_repunit, prop) %>% 
  unique() %>% 
  mutate(prop = replace_na(prop, 0),
         prop = round(prop, 2)) %>% 
  spread(repunit, prop) %>% 
  replace(is.na(.), 0)

write_csv(SA_rates, "./output/SelfAssignment_confusionMatrix.csv") #Output matrix

# ====

# Print corrected self-assignment rate for each reporting GMZ
sa_to_repu.maxLikelihood %>% 
  ungroup() %>% 
  group_by(repunit) %>% 
  summarize(correct = mean(repunit == inferred_repunit))

# Calculate mean scaled likelihood
sa_mean_scaled_like <- sa_to_repu %>% 
  mutate(repunit_f = factor(x = repunit, levels = RepGroups)) %>% #turn character into factor
  mutate(inferred_repunit_f = factor(x = inferred_repunit, levels = RepGroups)) %>% #turn character into factor
  group_by(repunit_f, inferred_repunit_f) %>% 
  summarise(mean_repu_scaled_like = mean(repu_scaled_like))  #take the mean of the scaled liklihoods

# Generate  heatmap  
ggplot(sa_mean_scaled_like, aes(x = repunit_f, y = inferred_repunit_f, z = mean_repu_scaled_like)) +
  geom_tile(aes(fill = mean_repu_scaled_like)) +
  geom_text(aes(label=sprintf(unlist(sa_mean_scaled_like[,3]), fmt='%#.2f'))) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), name = "Mean Scaled Likelihood\n") +
  xlab("Reporting Unit") +
  ylab("Inferred Reporting Unit") +
  theme_classic(base_size = 15)

ggsave("./output/Fig_selfAssignment.png", width = 7, height = 5, dpi = 300)


# =========================================

# 03 | 100% mixture simulations

# Create 100% pure mixture scenario for each GMZ
Proof100_scenarios_lwf <- lapply(RepGroups, function(x) tibble(repunit = x, ppn = 1.0))
# Add GMZs to labels
names(Proof100_scenarios_lwf) <- paste("All", RepGroups, sep = "-")

# Run 100% mixtures by leave-one-out analysis
Proof100_results_lwf <- assess_reference_loo(reference = refDat, 
                                             gen_start_col = 5, 
                                             reps = 50, # repetitions
                                             mixsize = 100, # number of samples
                                             alpha_repunit = Proof100_scenarios_lwf,
                                             alpha_collection = 10) # distribution of individuals from collections ~equal

# Analysis summary by GMZ level
Proof100_summary <- Proof100_results_lwf %>% 
  mutate(repunit_f = factor(x = repunit, levels = unlist(RepGroups))) %>%
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n)) %>% 
  filter(gsub(pattern="All-","",repunit_scenario) == repunit_f)

# Generate plot of pure 100% mixture simulation
Proof100_summary %>% 
  ggplot(aes(x = repunit_f, y = repprop_posterior_mean, fill = repunit_f)) +
  geom_boxplot() +
  scale_fill_manual(values = RepColors) +
  guides(fill = "none") +
  geom_hline(yintercept = 0.9, linetype = "dotted") +
  xlab("") +
  ylab("Posterior Mean Reporting Group Proportion") +
  theme_bw(base_size = 15)

ggsave("./output/Fig_100percent.png", width = 7, height = 5, dpi = 300)

# Calculate mean estimated GMZ proportion for each GMZ
Proof100_summary %>%
  group_by(repunit_f) %>%
  summarise(mean100 = mean(repprop_posterior_mean))


# =========================================

# 04 | Realistic fishery mixtures (variable composition)

# Set seed
set.seed(1)

# Define candidate GMZ mixture compositions at various alpha values
DirichletProps <- as.data.frame(rbind(
  MCMCpack::rdirichlet(1250, c(rep(0.01, 6))),
  MCMCpack::rdirichlet(1250, c(rep(0.1, 6))),
  MCMCpack::rdirichlet(1250, c(rep(1, 6))),
  MCMCpack::rdirichlet(1250, c(rep(10, 6)))))

# Define reporting group GMZs
colnames(DirichletProps) <- RepGroups

# Assign unique ID to each mixture candidate
DirichletProps$Analysis <- 1:nrow(DirichletProps)

# Create table with one GMZ reporting group per scenario
RFpropsDir <- DirichletProps %>%
  pivot_longer(-Analysis, names_to = "variable", values_to = "value")

# Keep four simulations per reporting group and 1% bin
RFkeep <- RFpropsDir %>%
  mutate(bin = pmin(floor(value * 100), 99)) %>%
  group_by(variable, bin) %>%
  slice_sample(n = 4)

# Keep all GMZ configurations for all scenario IDs
RFpropsDir <- RFpropsDir %>%
  filter(Analysis %in% unique(RFkeep$Analysis))

# Make a list of scenarios
RF_scenarios_lwf <- lapply(unique(RFpropsDir$Analysis), function(x) 
  RFpropsDir %>%
    filter(Analysis == x) %>%
    ungroup() %>%
    dplyr::select(-1) %>%
    setNames(., c("repunit", "ppn")))

# Assign unique ID name to each simulation scenario
names(RF_scenarios_lwf) <- paste("RF-", unique(RFpropsDir$Analysis), sep="")

# Leave-one-out analysis
RF_results_lwf <- assess_reference_loo(reference = refDat, 
                                       gen_start_col = 5, 
                                       reps = 1, 
                                       mixsize = 100,
                                       alpha_repunit = RF_scenarios_lwf,
                                       alpha_collection = 10,
                                       resampling_unit = "gene_copies",
                                       return_indiv_posteriors = F)

# Simulation output summary
RFsum_lwf <- RF_results_lwf %>%
  mutate(repunit_f = factor(x = repunit, levels = RepGroups)) %>%
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), reprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))

# Scatter plot comparing GMZ reporting group proportions
ggplot(RFsum_lwf, aes(x = repu_n_prop, y = reprop_posterior_mean, colour = repunit_f)) +
  geom_point() +
  scale_color_manual(values = RepColors) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0.1, slope = 1, lty=2) +
  geom_abline(intercept = -0.1, slope = 1, lty=2) +
  facet_wrap(~ repunit_f) +
  ylab("Posterior Mean Reporting Group Proportion") +
  guides(colour = "none") +
  xlab("True Mixture Proportion") +
  theme_bw(base_size = 15)

ggsave("./output/Fig_realisticMixtures.png", width = 8, height = 6, dpi = 300)

# How many estimates are within +/- 0.10 of the true proportion?
RFsum_lwf %>%
  group_by(repunit_f) %>%
  summarise(within0.1 = mean(abs(reprop_posterior_mean - repu_n_prop) <= 0.1),
            RMSE = sqrt(mean((reprop_posterior_mean - repu_n_prop)^2)))
