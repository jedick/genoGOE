# Load libraries and read data
library(ggplot2)
library(patchwork)

df <- read.csv("cleaned_data.csv", check.names = FALSE)

# Add a column for binary classification of Oxygen tolerance (Anaerobe or Other)
is_anaerobe <- grepl("anaerobe", df$`oxygen_tolerance.Oxygen tolerance`)
df$`Oxygen tolerance` <- ifelse(is_anaerobe, "Anaerobe", "Other")

# Create binary temperature groups based on median temperature
temp_median <- median(df$`culture_temp.Temperature`, na.rm = TRUE)
df$`Temperature group` <- ifelse(
  df$`culture_temp.Temperature` <= temp_median,
  "Low-T",
  "High-T"
)

# Density plot for GC content by Oxygen tolerance
p_oxygen <- ggplot(df, aes(x = `GC_content.GC-content`, fill = `Oxygen tolerance`)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(x = "GC content", y = "Density", fill = "Oxygen tolerance")

# Density plot for GC content by Temperature group
p_temp <- ggplot(df, aes(x = `GC_content.GC-content`, fill = `Temperature group`)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(x = "GC content", y = "Density", fill = "Temperature group")

# Arrange both plots in a single row
combined_plot <- p_oxygen + p_temp + plot_layout(nrow = 1)

ggsave("GC_content_visualization.png", combined_plot, width = 10, height = 4)
