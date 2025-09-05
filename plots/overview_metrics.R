# Working directory -------------------------------------------------------
setwd("E:/2024-07-UTBenchmark-SpC")

# Import libraries --------------------------------------------------------
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Read data ---------------------------------------------------------------
df <- read.csv("E:/2024-07-UTBenchmark-SpC/data/processed/plots/overview_metrics/Overview.csv", header = TRUE, sep=',') %>%
  filter(Sample != "Mouse_colon") %>%
  filter(Technology != "BMKMANU_S1000")

# Exclude 'Sample' and 'Technology' columns
cols_to_plot <- setdiff(names(df), c("Sample", "Technology"))

# Define colors
tech_colors <- c(
  "Visium" = "#f7c7ac",
  "Visium_HD_FFPE" = "#b5e6a2",
  "Visium_HD_FF" = "#fff07a",
  "BMKMANU_S1000" = "#a6c9ec",
  "BMKMANU_S3000" = "#e49edd",
  "Curio_Seeker" = "#d9d9d9"
)

# Function to create bar plots
make_bar_plot <- function(df, ycol) {
  # Create df_plot
  df_plot <- df %>%
    select(Sample, Technology, !!sym(ycol)) %>%
    rename(y_val = !!sym(ycol)) %>%
    mutate(
      label_text = case_when(
        is.na(y_val) ~ "NA",
        y_val > 1000 ~ scales::comma(y_val),
        TRUE ~ as.character(y_val)
      )
    )

  # Specify dodge width
  dodge_width <- 0.9
  
  # Determine if Percentage
  percentage = grepl("Percentage|Sequencing_saturation", ycol)
  
  # Determine y-axis formatting
  max_val <- max(df_plot$y_val, na.rm = TRUE)
  y_labels_func <- if (max_val > 1e6) {
    scales::label_number(scale = 1e-6, suffix = "M")
  } else {
    waiver()  # default formatting
  }

  # Plot
  p <- ggplot(df_plot, aes(x = Sample, y = y_val, fill = Technology)) +
    geom_bar(stat = "identity", position = position_dodge(width = dodge_width), color = "black") +
    geom_text(aes(label = label_text), position = position_dodge(width = dodge_width), hjust = -0.1, vjust = 0.5, size = 3.3) +
    coord_flip() +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(color = "black", size = 10),
      axis.title = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey30", linewidth = 0.25),
      axis.text = element_text(color = "black", size = 11),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    labs(title = ycol, y = NULL, x = NULL, fill = "Technology") +
    scale_fill_manual(values = tech_colors) +
    guides(fill = guide_legend(reverse = TRUE))
  
  p <- p + scale_y_continuous(
    breaks = if (percentage) seq(0, 101, by = 10) else pretty_breaks(n = 8),
    labels = y_labels_func,
    limits = if (percentage) c(0, 101) else NULL,
    expand = expansion(mult = c(0, 0.25))
  )
  
  height <- if (nrow(df) < 4) 1.5 else if (nrow(df) < 10) 3 else 4.5
  
  ggsave(paste0(gsub(" ", "_", ycol), ".png"), plot = p, width = 8, height = height)
}

# Loop over all columns
for (col_name in cols_to_plot) {
  make_bar_plot(df, col_name)
}
