library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)

load("../data/homoplasy_mutations.rda")
load("../data/homoplasy_mutations_norm.rda")

filtered_mutations <- homoplasy_nodes_annotated_byMutation %>%
  filter(synonym %in% c("phoR", "katG", "rpoB", "sigA"))

filtered_mutations_norm <- homoplasy_nodes_annotated_byMutation_norm %>%
  filter(synonym %in% c("phoR", "katG", "rpoB", "sigA"))

### RAW DATA ###

# Violin plot
raw_violin <- homoplasy_nodes_annotated_byMutation %>%
  ggplot(aes(x="", y=log2(RoHO))) +
  geom_violin(fill = "coral2") +
  geom_jitter(data = filtered_mutations, aes(color = synonym), alpha = 0.8, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.6) +
  scale_color_viridis(discrete = TRUE)
  labs(color = "Gene")
  
 # Histogram
raw_histogram <- homoplasy_nodes_annotated_byMutation %>%
  ggplot(aes(x = log2(RoHO))) +
  geom_histogram(binwidth = 0.25, fill="coral2", color="#e9ecef", alpha = 0.9)

# Boxplot
raw_boxplot <- homoplasy_nodes_annotated_byMutation %>%
  ggplot(aes(x="",y = log2(RoHO))) +
  geom_boxplot(fill = "coral2") +
  geom_jitter(data = filtered_mutations, aes(color = synonym), alpha = 0.8, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.6) +
  scale_color_viridis(discrete = TRUE) +
  labs(color = "Gene")


### NORM DATA ###

# Violin plot
norm_violin <- homoplasy_nodes_annotated_byMutation_norm %>%
  ggplot(aes(x="", y=log2(normalized_RoHO))) +
  geom_violin(fill = "coral2") +
  geom_jitter(data = filtered_mutations_norm, aes(color = synonym), alpha = 0.8, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.6) +
  scale_color_viridis(discrete = TRUE)
  labs(color = "Gene")

# Histogram
norm_histogram <- homoplasy_nodes_annotated_byMutation_norm %>%
  ggplot(aes(x = log2(normalized_RoHO))) +
  geom_histogram(binwidth = 0.25, fill="coral2", color="#e9ecef", alpha = 0.9)

# Boxplot
norm_boxplot <- homoplasy_nodes_annotated_byMutation_norm %>%
  ggplot(aes(x="",y = log2(normalized_RoHO))) +
  geom_boxplot(fill = "coral2") +
  geom_jitter(data = filtered_mutations_norm, aes(color = synonym), alpha = 0.8, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.6) +
  scale_color_viridis(discrete = TRUE) +
  labs(color = "Gene")


raw_violin + norm_violin
raw_histogram + norm_histogram
raw_boxplot + norm_boxplot  
  
