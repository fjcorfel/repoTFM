annotated_tree$accumulated_mutations_to_root <- 0
for (n_distance in seq_along(annotated_tree$distance_to_root)) {
  distance <- annotated_tree$distance_to_root[n_distance]
  
  filtered_nodes <- annotated_tree %>%
    filter(distance_to_root <= distance)
  
  accumulated_mutations <- sum(filtered_nodes$n_synonym_mutations)
  
  annotated_tree$accumulated_mutations_to_root[n_distance] <- accumulated_mutations
}



linear_model <- lm(accumulated_mutations_to_root ~ distance_to_root,
                   data = annotated_tree)

ggplot(annotated_tree, aes(x = distance_to_root, y = accumulated_mutations_to_root)) +
  geom_point(color = "red", size = 1.5, alpha = 0.6) +  # Puntos de datos
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +  # Línea de regresión
  labs(title = "N mutaciones sinónimas acumuladas vs Distancia a la Raíz",
       x = "Distancia a la Raíz",
       y = "Mutaciones Sinónimas Acumuladas") +
  ylim(-2000, 15000) +
  theme_minimal()

new_distances <- data.frame(distance_to_root = c(0.000206514))
predictions <- predict(linear_model, new_distances)
predictions
