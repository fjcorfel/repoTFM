library(ggplot2)
library(dplyr)

load("../data/homoplasy_nodes_genes.rda")

RoHO <- homoplasy_nodes$RoHO[homoplasy_nodes$RoHO < 20]
hist(RoHO, breaks = 1000)

homoplasy_nodes <- homoplasy_nodes %>%
  filter(RoHO < 20) %>%
  mutate(log2RoHO = log2(RoHO))


ggplot(homoplasy_nodes, aes(y = log2RoHO)) + 
  geom_boxplot(fill = "cyan", color = "black") +
  theme_minimal()


homoplasy_nodes_phoR <- homoplasy_nodes %>%
  filter(synonym == "phoR")

ggplot(homoplasy_nodes, aes(x ="", y=log2RoHO)) +
  geom_violin(fill = "lightgreen") +
  geom_jitter(data = homoplasy_nodes_phoR, color = "red")
