# Effect of probe count on gene-level log2 fold changes
#
# Usage:
# Rscript probe_count_log2fc.R [probe_set.csv] [output_root]

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

probe_set_path <- if (length(args) >= 1) args[[1]] else paste0(
  "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/",
  "data/processed/Visium_HD/SCA001/subsampled_100M/SpaceRanger/outs/probe_set.csv"
)
output_root <- if (length(args) >= 2) args[[2]] else paste0(
  "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/",
  "data/processed/plots/bam_plots/probe_count_effect"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

log2fc_paths <- c(
  "Mouse brain" = paste0(
    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/",
    "data/processed/plots/bam_plots/Mouse brain/gene_counts/",
    "Mouse brain_sample_size_100000000_pairwise_log2FC.csv"
  ),
  "Mouse liver" = paste0(
    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/",
    "data/processed/plots/bam_plots/Mouse liver/gene_counts/",
    "Mouse liver_sample_size_100000000_pairwise_log2FC.csv"
  )
)

comparison_technologies <- c(
  "Visium HD (FF)",
  "BMKMANU S3000",
  "Stereo-Seq",
  "Nova-ST",
  "Curio Seeker"
)

reference_analysis <- list(
  label = "Visium HD 3'",
  pattern = "^Visium HD 3'($| \\([0-9]+\\)$)",
  file_suffix = ""
)

probe_group_levels <- c("All genes", "1 probe", "2 probes", "3 probes", "4+ probes")
probe_group_colors <- c(
  "All genes" = "#999999",
  "1 probe" = "#0072B2",
  "2 probes" = "#E69F00",
  "3 probes" = "#009E73",
  "4+ probes" = "#CC79A7"
)
probe_group_legend_labels <- c(
  "All genes" = "All genes",
  "1 probe" = "Genes targeted by 1 probe",
  "2 probes" = "Genes targeted by 2 probes",
  "3 probes" = "Genes targeted by 3 probes",
  "4+ probes" = "Genes targeted by 4+ probes"
)

assign_probe_group <- function(number_of_probes) {
  case_when(
    number_of_probes == 1 ~ "1 probe",
    number_of_probes == 2 ~ "2 probes",
    number_of_probes == 3 ~ "3 probes",
    number_of_probes >= 4 ~ "4+ probes",
    TRUE ~ NA_character_
  )
}

# The metadata lines begin with '#'. The gene symbol is the middle component
# of probe_id: gene_id|gene_symbol|probe_hash.
included_probes <- read_csv(
  probe_set_path,
  comment = "#",
  show_col_types = FALSE
) %>%
  filter(included) %>%
  mutate(gene = str_split_fixed(probe_id, fixed("|"), 3)[, 2])

probe_counts <- included_probes %>%
  count(gene_id, gene, name = "number_of_probes") %>%
  mutate(
    probe_group = assign_probe_group(number_of_probes),
    probe_group = factor(
      probe_group,
      levels = probe_group_levels[-1]
    )
  )

probe_group_summary <- probe_counts %>%
  count(probe_group, name = "number_of_genes", .drop = FALSE) %>%
  mutate(
    percentage_of_included_genes = 100 * number_of_genes / sum(number_of_genes)
  )

write_csv(probe_counts, file.path(output_root, "included_gene_probe_counts.csv"))
write_csv(probe_group_summary, file.path(output_root, "probe_count_group_summary.csv"))

probe_inventory_plot <- ggplot(
  probe_group_summary,
  aes(x = probe_group, y = number_of_genes, fill = probe_group)
) +
  geom_col(width = 0.72, color = "black", linewidth = 0.4) +
  geom_text(
    aes(
      label = format(
        number_of_genes,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
      )
    ),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_fill_manual(
    values = probe_group_colors[probe_group_levels[-1]],
    guide = "none",
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = function(x) format(
      x,
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    x = "Included probes per gene",
    y = "Number of genes"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(linewidth = 0.4, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey70", linewidth = 0.25),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12)
  )

ggsave(
  file.path(output_root, "probe_count_group_gene_counts.png"),
  plot = probe_inventory_plot,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)

# The log2FC table is keyed by gene symbol. Exclude symbols assigned to more
# than one gene ID because their number of probes is not uniquely defined.
ambiguous_symbols <- probe_counts %>%
  count(gene, name = "number_of_gene_ids") %>%
  filter(number_of_gene_ids > 1) %>%
  left_join(probe_counts, by = "gene") %>%
  arrange(gene, gene_id)

write_csv(
  ambiguous_symbols,
  file.path(output_root, "ambiguous_probe_gene_symbols.csv")
)

unambiguous_probe_counts <- probe_counts %>%
  anti_join(ambiguous_symbols %>% distinct(gene), by = "gene")

analyze_reference <- function(tissue, pairwise_log2fc, reference_analysis, output_dir) {
  reference_label <- reference_analysis$label
  reference_pattern <- reference_analysis$pattern
  file_suffix <- reference_analysis$file_suffix

  # Pairwise comparisons occur in only one direction in the input table.
  # Reorient them so log2FC is always comparison/reference: retain the sign
  # when the reference is tech2 and reverse it when the reference is tech1.
  oriented_log2fc <- bind_rows(
    pairwise_log2fc %>%
      filter(
        str_detect(tech2, reference_pattern),
        !str_detect(tech1, reference_pattern)
      ) %>%
      transmute(
        comparison_technology = tech1,
        reference_technology = tech2,
        pair_id = paste(tech1, tech2, sep = " versus "),
        gene,
        log2FC
      ),
    pairwise_log2fc %>%
      filter(
        str_detect(tech1, reference_pattern),
        !str_detect(tech2, reference_pattern)
      ) %>%
      transmute(
        comparison_technology = tech2,
        reference_technology = tech1,
        pair_id = paste(tech1, tech2, sep = " versus "),
        gene,
        log2FC = -log2FC
      )
  ) %>%
    filter(!is.na(log2FC)) %>%
    mutate(
      # Combine replicated datasets such as Visium HD 3' (1) and (2) under
      # one technology label while retaining their distinct pair IDs.
      technology = str_remove(comparison_technology, " \\([0-9]+\\)$")
    ) %>%
    filter(technology %in% comparison_technologies)

  if (nrow(oriented_log2fc) == 0) {
    stop("No comparisons against ", reference_label, " found for ", tissue)
  }

  expected_pair_counts <- oriented_log2fc %>%
    distinct(technology, pair_id) %>%
    count(technology, name = "expected_number_of_pairwise_comparisons")

  # Average replicate-level comparisons and require each retained gene to be
  # present in every pairwise comparison available for that technology.
  summarized_log2fc <- oriented_log2fc %>%
    group_by(technology, gene) %>%
    summarise(
      log2FC = mean(log2FC),
      number_of_pairwise_comparisons = n_distinct(pair_id),
      .groups = "drop"
    ) %>%
    left_join(expected_pair_counts, by = "technology") %>%
    filter(
      number_of_pairwise_comparisons == expected_number_of_pairwise_comparisons
    ) %>%
    mutate(
      tissue = tissue,
      reference_technology = reference_label,
      technology = factor(technology, levels = comparison_technologies)
    )

  analysis_data <- summarized_log2fc %>%
    inner_join(
      unambiguous_probe_counts %>%
        select(gene, gene_id, number_of_probes, probe_group),
      by = "gene"
    )

  comparison_group_counts <- analysis_data %>%
    count(technology, probe_group, name = "number_of_genes", .drop = FALSE)

  write_csv(
    analysis_data,
    file.path(
      output_dir,
      paste0("probe_count_pairwise_log2fc", file_suffix, ".csv")
    )
  )
  write_csv(
    comparison_group_counts,
    file.path(
      output_dir,
      paste0("probe_count_log2fc_group_counts", file_suffix, ".csv")
    )
  )

  # Duplicate each comparison once to add its complete gene distribution next
  # to the four probe-count strata. Values remain on the original log2FC scale.
  plot_data <- bind_rows(
    analysis_data %>% mutate(plot_group = "All genes"),
    analysis_data %>% mutate(plot_group = as.character(probe_group))
  ) %>%
    mutate(plot_group = factor(plot_group, levels = probe_group_levels))

  plot <- ggplot(
    plot_data,
    aes(x = technology, y = log2FC, fill = plot_group)
  ) +
    geom_vline(
      xintercept = 1.5,
      color = "grey50",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
    geom_boxplot(
      position = position_dodge2(width = 0.85, preserve = "single"),
      width = 0.8,
      color = "black",
      linewidth = 0.4,
      outlier.shape = NA
    ) +
    scale_fill_manual(
      values = probe_group_colors,
      limits = probe_group_levels,
      labels = probe_group_legend_labels[probe_group_levels],
      drop = FALSE
    ) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(add = c(0.5, 0.5))
    ) +
    scale_y_continuous(
      breaks = seq(-6, 6, by = 1),
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(ylim = c(-6, 6)) +
    labs(
      x = NULL,
      y = "Log2FC relative to Visium HD 3'",
      fill = "Gene set"
    ) +
    theme_minimal() +
    theme(
      axis.line.x = element_line(linewidth = 0.4, color = "black"),
      axis.line.y = element_line(linewidth = 0.4, color = "black"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey70", linewidth = 0.25),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(color = "black", size = 12),
      axis.ticks.y = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length.y = grid::unit(2, "mm"),
      legend.position = "right",
      plot.margin = margin(10, 20, 10, 10)
    )

  ggsave(
    file.path(
      output_dir,
      paste0("probe_count_log2fc_boxplots", file_suffix, ".png")
    ),
    plot = plot,
    width = 11,
    height = 6,
    dpi = 300,
    bg = "white"
  )

  message(
    tissue,
    ": reference comparison(s): ",
    reference_label
  )
  print(comparison_group_counts)
}

analyze_tissue <- function(tissue, log2fc_path) {
  output_dir <- file.path(output_root, tissue)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pairwise_log2fc <- read_csv(log2fc_path, show_col_types = FALSE)
  analyze_reference(
    tissue,
    pairwise_log2fc,
    reference_analysis,
    output_dir
  )
  message("Saved results to ", output_dir)
}

print(probe_group_summary)
for (tissue in names(log2fc_paths)) {
  analyze_tissue(tissue, log2fc_paths[[tissue]])
}
