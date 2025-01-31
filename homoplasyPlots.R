library(ggplot2)
library(dplyr)

load("../data/homoplasy_mutations.rda")

# Violin plot
phor_mutations <- homoplasy_nodes_annotated_byMutation %>%
  filter(synonym == "phoR")

homoplasy_nodes_annotated_byMutation %>%
  ggplot(aes(x="", y=log2(RoHO))) +
  geom_violin(fill = "coral2") +
  geom_point(data = phor_mutations, color = "darkslateblue", alpha = 0.8, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.6)
  
 # Histogram
homoplasy_nodes_annotated_byMutation %>%
  ggplot(aes(x = log2(RoHO))) +
  geom_histogram(binwidth = 0.1, fill="coral2", color="#e9ecef", alpha = 0.9)

# Boxplot
homoplasy_nodes_annotated_byMutation %>%
  ggplot(aes(x="",y = log2(RoHO))) +
  geom_boxplot(fill = "coral2") +
  geom_jitter(data = phor_mutations, color = "darkslateblue", alpha = 0.8, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.6)
  
