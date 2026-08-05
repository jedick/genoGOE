#!/usr/bin/env Rscript

# Load libraries and read data
library(ggplot2)
library(patchwork)

df <- read.csv("cleaned_data_with_Zc.csv", check.names = FALSE)

# Binary classification of Oxygen tolerance (Anaerobe or Other)
is_anaerobe <- grepl("anaerobe", df$`oxygen_tolerance.Oxygen tolerance`, ignore.case = TRUE)
df$`Oxygen tolerance` <- ifelse(is_anaerobe, "Anaerobe", "Other")

# Keep only rows with all variables needed for plotting
plot_df <- subset(
  df,
  !is.na(`culture_temp.Temperature`) &
    !is.na(`GC_content.GC-content`) &
    !is.na(Zc)
)

# Scatter: GC content vs temperature with linear fits by Oxygen tolerance
p_gc <- ggplot(
  plot_df,
  aes(
    x = `culture_temp.Temperature`,
    y = `GC_content.GC-content`,
    color = `Oxygen tolerance`
  )
) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  theme_minimal() +
  labs(
    x = "Culture temperature (°C)",
    y = "GC content (%)",
    color = "Oxygen tolerance",
    title = "A"
  ) +
  scale_color_manual(values = c("Anaerobe" = "2", "Other" = "4")) +
  theme(plot.title = element_text(face = "bold"))

# Scatter: Zc vs temperature with linear fits by Oxygen tolerance
p_zc <- ggplot(
  plot_df,
  aes(
    x = `culture_temp.Temperature`,
    y = Zc,
    color = `Oxygen tolerance`
  )
) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  theme_minimal() +
  labs(
    x = "Culture temperature (°C)",
    y = "Zc",
    color = "Oxygen tolerance",
    title = "B"
  ) +
  scale_color_manual(values = c("Anaerobe" = "2", "Other" = "4")) +
  theme(plot.title = element_text(face = "bold"))

# Arrange plots in a single row with shared legend
combined_plot <- p_gc + p_zc + plot_layout(nrow = 1, guides = "collect")

combined_plot <- combined_plot & theme(legend.position = "bottom")

ggsave("GC_Zc_vs_temperature.png", combined_plot, width = 10, height = 4, dpi = 300)

