# =========================================
# Lake whitefish GSI figures: STRUCTURE, RUBIAS mixture proportions, and PCA
# =========================================

library(tidyverse)
library(adegenet)

RepColors <- c("#A9D18E", "darkgreen", "#C48AC1", "#9DC3E6", "#DBAF5C", "#ADB9CA")

BSMod <- read_csv("./output/GSI_estimates_BSMod_SeasonYearZone.csv")
BSMod_year <- read_csv("./output/GSI_estimates_BSMod_YearZone.csv")

# =========================================

# 01 | STRUCTURE (Figure 2): Ancestry proportions for west-side Green Bay area, K = 2 
# NOTE: Q1 = lineage L1, Q2 = lineage L2

# Read Q values
q2 <- read_csv("./output/NW_K2_MeanQ.csv") %>%
  select(SampleID = Sample_ID, cluster = variable_new, prob = MeanQ, Location = Collection) %>%
  mutate(Location = as.factor(Location)) 

# Plot STRUCTURE clustering 
q2 %>%
  ggplot(aes(SampleID, prob, fill = factor(cluster))) +
  geom_col(color = "gray", width = 1) +
  facet_grid(~fct_relevel(Location, "Menominee River", "Riley Bay"), 
             switch = "x", 
             scales = "free", 
             space = "free") +
  labs(x = "Cluster", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_fill_manual(values = c("#E1BE6A", "#40B0A6")) +
  theme_minimal(base_size = 16) + 
  theme(panel.spacing.x = unit(0.05, "lines"),
        axis.text.x = element_blank(),
        panel.grid = element_blank()) + 
  guides(fill=guide_legend(title = "Genetic lineage")) +
  theme(strip.text.x = element_text(angle = -60)) 

ggsave("./output/Fig2_STRUCTURE.png", width = 12, height = 4, dpi = 300) # Save plot


# =========================================

# 02 | Seasonal mixture proportions (Figure 4): Years 2020 and 2022

BSMod %>%
  separate(mixture_collection, c("MZ", "Year", "Season"), sep = " ") %>%
  filter(str_detect(Year, "2020|2022")) %>%
  filter(str_detect(MZ, "WFM-00|WFM-01|WFM-02|WFM-06|WFM-08|WI-01|WI-02|WI-03")) %>%
  mutate(Season = fct_relevel(Season, "Spring", "Summer", "Fall")) %>%
  ggplot(aes(x = Season, y = bs_corrected_repunit_ppn, fill = fct(repunit))) +
  geom_col(color = "black") +
  facet_grid(Year ~ MZ) +
  labs(y="Boostrap-corrected mixture proportions") +
  labs(fill = "Reporting Unit") +
  scale_fill_manual(values = RepColors) +
  scale_y_continuous(expand = c(0,.01)) +
  theme_bw() +
  theme(text = element_text(size = 13),
        axis.text.x = element_text(angle = 25, margin = margin(t = 6)))

ggsave("./output/Fig4_seasonal_GSI.png", width = 12, height = 5, dpi = 300) # Save plot


# =========================================

# 03 | Long term mixture proportions (Figure 5)

BSMod_year %>%
  separate(mixture_collection, c("MZ", "Year"), sep = " ") %>%
  filter(str_detect(MZ, "WI-01|WI-02|WI-03|WFM-08")) %>%
  ggplot(aes(x = Year, y = bs_corrected_repunit_ppn, fill = repunit)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  labs(y="Boostrap-corrected mixture proportions") +
  labs(fill = "Reporting Unit") +
  facet_wrap(~MZ, scales = "free") +
  scale_fill_manual(values = RepColors) +
  scale_y_continuous(expand = c(0,.01)) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 25, margin = margin(t = 6)))

ggsave("./output/Fig5_longterm_GSI.png", width = 7, height = 7, dpi = 300) # Save plot


# =========================================

# 04 | Extra management zones mixture proportions (Suppl Fig 4)

BSMod_year %>%
  separate(mixture_collection, c("MZ", "Year"), sep = " ") %>%
  filter(str_detect(MZ, "WFM-03|WFM-04|WFM-05")) %>%
  ggplot(aes(x = Year, y = bs_corrected_repunit_ppn, fill = repunit)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  labs(y="Boostrap-corrected mixture proportions") +
  labs(fill = "Reporting Unit") +
  facet_wrap(~MZ, scales = "free") +
  scale_fill_manual(values = RepColors) +
  scale_y_continuous(expand = c(0,.01)) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, margin = margin(t = 6)))

ggsave("./output/Fig5_Suppl_GSI.png", width = 7, height = 7, dpi = 300) # Save plot




