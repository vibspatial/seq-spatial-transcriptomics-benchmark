# Gene-biotype overview plots
# Run: Rscript gene_biotype_plots.R [directory containing the overview CSVs]

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else paste0(
  "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/",
  "data/processed/plots/genes/biotypes"
)

biotype_summary <- read_csv(
  file.path(data_dir, "gene_biotype_summary_all_datasets.csv"),
  show_col_types = FALSE
)
if (!"biotype_group" %in% names(biotype_summary)) {
  stop(
    "The summary table has no 'biotype_group' column. Rerun ",
    "gene_biotype_summary.ipynb before creating the figures."
  )
}

# Exact ordering used by overview_metrics.R. With coord_flip(), the first
# sample type is displayed at the bottom and the last at the top.
technology_levels <- c(
  "Visium HD 3'", "Visium HD (FF)", "Visium HD (FFPE)", "Visium (FFPE)",
  "Stereo-Seq", "Nova-ST", "Curio Seeker", "BMKMANU S3000", "BMKMANU S1000"
)
sample_type_levels <- c(
  "Rice root", "Human lung", "Mouse lung", "Mouse liver", "Mouse brain"
)

biotype_colors <- c(
  "Protein-coding" = "#0072B2",
  "IG genes" = "#CC79A7",
  "TR genes" = "#D55E00",
  "lncRNA" = "#009E73",
  "rRNA" = "#E69F00",
  "tRNA" = "#56B4E9",
  "Pseudogene" = "#F0E442",
  "Other ncRNA" = "#332288"
)
biotype_order <- c(
  "Protein-coding", "IG genes", "TR genes", "lncRNA", "rRNA", "tRNA",
  "Pseudogene", "Other ncRNA"
)

format_replicate <- function(x) {
  ifelse(is.na(x), "", paste0(" (", format(x, trim = TRUE, scientific = FALSE), ")"))
}

prepare_metadata <- function(df) {
  df %>% mutate(
    technology = factor(technology, levels = technology_levels),
    sample_type = factor(sample_type, levels = sample_type_levels),
    dataset_label = paste0(as.character(technology), format_replicate(replicate))
  )
}

# Fixed positions recreate the technology dodge within every sample type while
# still allowing each gene-biotype bar to be stacked.
make_dataset_positions <- function(df) {
  datasets <- df %>%
    prepare_metadata() %>%
    distinct(sample_type, technology, replicate, dataset_label) %>%
    arrange(sample_type, technology, coalesce(replicate, 0))

  maximum_datasets_per_sample <- datasets %>%
    count(sample_type) %>%
    summarise(maximum = max(n)) %>%
    pull(maximum)
  bar_width <- 0.9 / maximum_datasets_per_sample

  datasets %>%
    group_by(sample_type) %>%
    mutate(
      dataset_number = row_number(), number_of_datasets = n(),
      sample_number = as.numeric(sample_type),
      bar_position = sample_number +
        (dataset_number - (number_of_datasets + 1) / 2) * bar_width
    ) %>%
    ungroup() %>%
    mutate(bar_width = .env$bar_width)
}

dataset_positions <- make_dataset_positions(biotype_summary)
shade_df <- tibble(
  sample_type = factor(sample_type_levels, levels = sample_type_levels),
  sample_number = seq_along(sample_type_levels)
) %>% mutate(grey_bg = sample_type %in% c("Mouse brain", "Mouse lung", "Rice root"))

# Styling copied from make_bar_plot() in overview_metrics.R.
overview_theme <- function() {
  theme_minimal() + theme(
    axis.line.x = element_line(linewidth = 0.4),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(color = "black", size = 10),
    axis.title.x = element_text(size = 12), axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey30", linewidth = 0.25),
    axis.text = element_text(color = "black", size = 11),
    legend.position = "right", plot.margin = margin(10, 20, 10, 10)
  )
}

add_sample_type_backgrounds <- function(plot) {
  plot + geom_rect(
    data = shade_df,
    aes(xmin = sample_number - 0.5, xmax = sample_number + 0.5,
        ymin = -Inf, ymax = Inf, alpha = grey_bg),
    inherit.aes = FALSE, fill = "grey35", color = NA
  ) + scale_alpha_manual(values = c(`TRUE` = 0.12, `FALSE` = 0), guide = "none")
}

add_sample_type_group_lines <- function(plot) {
  plot + geom_segment(
    data = shade_df,
    aes(
      x = sample_number - 0.45,
      xend = sample_number + 0.45,
      y = -30,
      yend = -30
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.4
  )
}

add_common_scales <- function(plot, x_title) {
  plot + coord_flip(clip = "off") +
    scale_x_continuous(
      breaks = seq_along(sample_type_levels), labels = sample_type_levels,
      limits = c(0.5, length(sample_type_levels) + 0.5),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 10),
      minor_breaks = seq(5, 95, by = 10),
      limits = c(-35, 100),
      expand = expansion(mult = c(0, 0.1))
    ) +
    labs(x = NULL, y = str_wrap(x_title, width = 35)) + overview_theme()
}

save_plot <- function(plot, filename) {
  ggsave(
    file.path(data_dir, paste0(filename, ".png")), plot = plot,
    width = 9.5, height = 7, dpi = 300, bg = "white"
  )
}

make_biotype_plot <- function(metric, x_title, filename) {
  plot_data <- biotype_summary %>%
    prepare_metadata() %>%
    filter(biotype_group %in% biotype_order) %>%
    mutate(
      biotype_group = factor(biotype_group, levels = biotype_order)
    ) %>%
    group_by(sample_type, technology, replicate, dataset_label, biotype_group) %>%
    summarise(value = sum(.data[[metric]], na.rm = TRUE), .groups = "drop") %>%
    left_join(
      dataset_positions,
      by = c("sample_type", "technology", "replicate", "dataset_label")
    )

  plot <- ggplot()
  plot <- add_sample_type_backgrounds(plot)
  plot <- add_sample_type_group_lines(plot)
  plot <- plot +
    geom_col(
      data = plot_data,
      aes(x = bar_position, y = value, fill = biotype_group),
      width = unique(plot_data$bar_width), color = NA,
      position = position_stack(reverse = TRUE)
    ) +
    geom_rect(
      data = dataset_positions,
      aes(
        xmin = bar_position - bar_width / 2,
        xmax = bar_position + bar_width / 2,
        ymin = 0,
        ymax = 100
      ),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.5
    ) +
    geom_text(
      data = dataset_positions,
      aes(x = bar_position, y = -1.5, label = dataset_label),
      hjust = 1, vjust = 0.5, size = 2.6, color = "black"
    ) +
    scale_fill_manual(
      values = biotype_colors,
      limits = biotype_order,
      drop = FALSE
    ) +
    labs(fill = "Gene biotype")

  save_plot(add_common_scales(plot, x_title), filename)
}

make_biotype_plot(
  "percentage_of_counts", "Contribution to unique UMI counts in tissue (%)",
  "gene_biotype_count_contribution_overview"
)
make_biotype_plot(
  "percentage_of_detected_features", "Percentage of genes detected in tissue (%)",
  "gene_biotype_detected_features_overview"
)
