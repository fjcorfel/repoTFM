library(ggplot2)
library(dplyr)
library(viridis)

load("../data/homoplasy_mutations.rda")

filtered_mutations <- homoplasy_nodes_annotated_byMutation %>%
  filter(synonym %in% c("phoR", "katG", "rpoB", "sigA"))


# Violin plot
homoplasy_nodes_annotated_byMutation %>%
  ggplot(aes(x="", y=log2(RoHO))) +
  geom_violin(fill = "coral2") +
  geom_jitter(data = filtered_mutations, aes(color = synonym), alpha = 0.8, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.6) +
  scale_color_viridis(discrete = TRUE)
  labs(color = "Gene")
  
 # Histogram
homoplasy_nodes_annotated_byMutation %>%
  ggplot(aes(x = log2(RoHO))) +
  geom_histogram(binwidth = 0.25, fill="coral2", color="#e9ecef", alpha = 0.9)

# Boxplot
homoplasy_nodes_annotated_byMutation %>%
  ggplot(aes(x="",y = log2(RoHO))) +
  geom_boxplot(fill = "coral2") +
  geom_jitter(data = filtered_mutations, aes(color = synonym), alpha = 0.8, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.6) +
  scale_color_viridis(discrete = TRUE) +
  labs(color = "Gene")
  
  
  
