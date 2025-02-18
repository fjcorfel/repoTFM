library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)
library(plotly)
library(hrbrthemes)


load("../data/homoplasy_mutations.rda")

# Histogram with RoHO distribution
histogram <- homoplasy_mutations %>%
  ggplot(aes(x = log2(RoHO))) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.4) +
  ggtitle("Distribution of RoHO values")

histogram

#------------------------------------------------------------------------------

# RoHO distribution with control genes
violin <- homoplasy_mutations %>%
  ggplot(aes(x="", y=log2(RoHO))) +
  geom_violin(fill = "coral2", color = "#e9ecef",trim = FALSE, alpha = 0.8) +
  geom_jitter(data = homoplasy_mutations %>%
                filter(synonym %in% c("phoR", "katG", "rpoB", "sigA")),
              aes(color = synonym), alpha = 0.8, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.4) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_color_viridis(discrete = TRUE)

boxplot <- homoplasy_mutations %>%
  ggplot(aes(x="", y=log2(RoHO))) +
  geom_boxplot(fill = "coral2", alpha = 0.8) +
  geom_jitter(data = homoplasy_mutations %>%
                filter(synonym %in% c("phoR", "katG", "rpoB", "sigA")),
              aes(color = synonym), alpha = 0.8, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.4) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_viridis(discrete = TRUE)

violin + boxplot + plot_annotation(title = "Distribution of RoHO values")
