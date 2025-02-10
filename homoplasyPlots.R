library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)
library(plotly)
library(hrbrthemes)


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
load("../data/top_RoHO_mutations_v2.0.rda")
load("../data/homoplasy_mutations.rda")
load("../data/n_nodes_per_homoplasy.rda")

data <- homoplasy_nodes_annotated_byMutation %>%
  left_join(n_homoplasy_nodes, by = "mutation")

scatter_nodes_RoHO_all <- data %>%
  ggplot(aes(x = n_nodes, y = RoHO)) +
  geom_point(color = "#69b3a2", alpha = 0.8, size = 1) +
  geom_point(data = top_RoHO_mutations %>% filter(pvalue <= 0.05), color = "#B3697A", alpha = 0.8,size = 1) +
  ggtitle("Relation between number of nodes per mutation and RoHO value") + 
  xlab("Number of nodes per mutation") +
  ylab("RoHO")

scatter_nodes_RoHO_top <- top_RoHO_mutations %>%
  ggplot(aes(x = n_nodes, y = RoHO)) +
  geom_point(color = "#69b3a2", alpha = 0.8, size = 1) +
  geom_point(data = top_RoHO_mutations %>% filter(pvalue <= 0.05),  color = "#B3697A", alpha = 0.8,size = 1) +
  ggtitle("Relation between number of nodes per mutation and RoHO value") + 
  xlab("Number of nodes per mutation") +
  ylab("RoHO")


ggplotly(scatter_nodes_RoHO_all, tooltip = "text")
ggplotly(scatter_nodes_RoHO_top, tooltip = "text")
