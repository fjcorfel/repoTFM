library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)

load("../data/homoplasy_mutations.rda")
load("../data/homoplasy_mutations_norm.rda")

raw_data <- homoplasy_nodes_annotated_byMutation %>%
  mutate(Type = "Raw") %>%
  select(mutation,
         synonym,
         Rv_number,
         RoHO,
         Type)

norm_data <- homoplasy_nodes_annotated_byMutation_norm %>%
  mutate(Type = "Normalized") %>%
  select(mutation,
         synonym,
         Rv_number,
         norm_RoHO,
         Type) %>%
  rename(RoHO = norm_RoHO)

RoHO_data <- bind_rows(raw_data, norm_data)

filtered_mutations <- RoHO_data %>%
  filter(synonym %in% c("phoR", "katG", "rpoB", "sigA"))


### RAW DATA ###

# Violin plot
violin <- RoHO_data %>%
  ggplot(aes(x=Type, y=log2(RoHO), fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Use full dataset
  geom_jitter(data = filtered_mutations, aes(color = synonym), alpha = 0.8, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.6) +
  labs(color = "Gene") +
  scale_color_viridis(discrete = TRUE)

  
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
  labs(color = "Gene") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())



# N Nodes Vs RoHO
load("../data/top_RoHO_mutations.rda")
load("../data/homoplasy_mutations.rda")
load("../data/n_nodes_per_homoplasy.rda")

data <- homoplasy_nodes_annotated_byMutation %>%
  left_join(n_homoplasy_nodes, by = "mutation")

data %>%
  ggplot(aes(x = n_nodes, y = RoHO)) +
  geom_point(alpha = 0.8, size = 1) +
  geom_point(data = top_RoHO_mutations, aes(x = n_nodes, y = RoHO, color = pvalue <= 0.05)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  geom_hline(yintercept = min(top_RoHO_mutations$RoHO),
             linetype = "dashed",
             color = "blue",
             alpha = 0.5) +
  coord_cartesian(ylim = c(0, 100))
