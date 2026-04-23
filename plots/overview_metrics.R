# Set
setwd("/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/plots/overview_metrics")

# Import libraries --------------------------------------------------------
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Read data ---------------------------------------------------------------
df <- read.csv("/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/plots/overview_metrics/Overview.csv", header = TRUE, sep=',') %>%
  filter(Sample != "Mouse_colon") %>%
  mutate(Technology = recode(
          Technology,
          "Visium_HD_3" = "Visium HD 3'",
          "Visium_HD_FF" = "Visium HD (FF)",
          "Visium_HD_FFPE" = "Visium HD (FFPE)",
          "Visium" = "Visium (FFPE)",
          "Stereo-Seq" = "Stereo-Seq",
          "Nova-ST" = "Nova-ST",
          "Curio_Seeker" = "Curio Seeker",
          "BMKMANU_S3000" = "BMKMANU S3000",
          "BMKMANU_S1000" = "BMKMANU S1000"),
        Technology = factor(
          Technology, levels = c(
          "Visium HD 3'",
          "Visium HD (FF)",
          "Visium HD (FFPE)",
          "Visium (FFPE)",
          "Stereo-Seq",
          "Nova-ST",
          "Curio Seeker",
          "BMKMANU S3000",
          "BMKMANU S1000")),
        Sample = recode(Sample,
          "Rice_root" = "Rice root",
          "Human_lung" = "Human lung",
          "Mouse_lung" = "Mouse lung",
          "Mouse_liver" = "Mouse liver",
          "Mouse_brain" = "Mouse brain"),
        Sample = factor(Sample, levels = c(
          "Rice root",
          "Human lung",
          "Mouse lung",
          "Mouse liver",
          "Mouse brain"
          )))

# Exclude 'Sample', 'Technology' and 'Replicate' columns
cols_to_plot <- setdiff(names(df), c("Sample", "Technology", "Replicate"))

# Specify pretty names
pretty_names <- c(
  Spot_diameter_um = "Spot diameter (µm)",
  Total_capture_area_mm2_reported = "Capture area (mm²) - reported",
  Total_capture_area_mm2_measured = "Capture area (mm²) - measured",
  
  Genes_targeted = "Genes targeted",
  
  Total_number_of_spots = "Spots",
  Number_of_spots_in_tissue = "Spots located in tissue",
  Percentage_of_spots_in_tissue = "Percentage of spots located in tissue (%)",
  
  Capture_area_covered_by_tissue_mm2 = "Capture area covered by tissue (mm²)",
  Percentage_capture_area_covered_by_tissue = "Percentage of capture area covered by tissue (%)",
  Percentage_of_tissue_covered_by_spots = "Percentage of tissue covered by spots (%)",
  Effective_capture_area_in_tissue_mm2 = "Effective capture area in tissue (mm²)",
  
  Q30_bases_in_spatial_barcode = "Q30 bases in spatial barcode",
  Q30_bases_in_probe_read = "Q30 bases in probe/read",
  Q30_bases_in_UMI = "Q30 bases in UMI",
  
  Number_of_spots_with_at_least_1_UMI_in_tissue = "Spots in tissue with ≥1 UMI count",
  Percentage_of_spots_recovered_in_tissue = "Percentage of spots recovered in tissue (%)",
  
  Total_number_of_reads = "Reads",
  Total_number_of_genes = "Genes",
  Total_number_of_UMIs = "UMI counts",
  Total_number_of_unique_UMIs = "Unique UMI counts",
  
  Number_of_genes_in_tissue = "Genes detected in tissue",
  Number_of_unique_UMIs_in_tissue = "Unique UMI counts detected in tissue",
  Percentage_of_genes_in_tissue = "Percentage of genes detected in tissue (%)",
  Percentage_of_targeted_genes_found_in_capture_area = "Percentage of targeted genes detected in capture area (%)",
  Percentage_of_targeted_genes_found_in_tissue = "Percentage of targeted genes detected in tissue (%)",
  Percentage_of_unique_UMIs_in_tissue = "Percentage of UMI counts detected in tissue (%)",
  
  Sequencing_saturation_reported = "Sequencing saturation (%) - reported",
  Sequencing_saturation_calculated_total_reads = "Sequencing saturation (%)",
  Sequencing_saturation_calculated_total_UMIs = "Sequencing saturation (%)",
  
  Unique_UMIs_per_mm2_in_tissue = "Unique UMI counts detected in tissue/mm²",
  
  Percentage_of_reads_that_contains_a_usable_UMI = "Percentage of reads contributing to UMI counts (%)",
  ratio_reads_unique_UMIs = "Reads/unique UMI",
  
  Number_of_spots_with_at_least_1_UMI_in_tissue_100M = "Spots in tissue with ≥1 UMI count",
  Percentage_of_spots_recovered_in_tissue_100M = "Percentage of spots recovered in tissue (%)",
  Total_number_of_reads_100M = "Reads",
  Total_number_of_genes_100M = "Genes",
  Total_number_of_UMIs_100M = "UMI counts",
  Total_number_of_unique_UMIs_100M = "Unique UMI counts",
  Number_of_genes_in_tissue_100M = "Genes detected in tissue",
  Number_of_unique_UMIs_in_tissue_100M = "Unique UMI counts detected in tissue",
  Percentage_of_genes_in_tissue_100M = "Percentage of genes detected in tissue (%)",
  Percentage_of_targeted_genes_found_in_capture_area_100M = "Percentage of targeted genes detected in capture area (%)",
  Percentage_of_targeted_genes_found_in_tissue_100M = "Percentage of targeted genes detected in tissue (%)",
  Percentage_of_unique_UMIs_in_tissue_100M = "Percentage of UMI counts detected in tissue (%)",
  Sequencing_saturation_reported_100M = "Sequencing saturation (%) - reported",
  Sequencing_saturation_calculated_total_reads_100M = "Sequencing saturation (%)",
  Sequencing_saturation_calculated_total_UMIs_100M = "Sequencing saturation (%)",
  Unique_UMIs_per_mm2_in_tissue_100M = "Unique UMI counts detected in tissue/mm²",
  Percentage_of_reads_that_contains_a_usable_UMI_100M = "Percentage of reads contributing to UMI counts (%)",
  ratio_reads_unique_UMIs_100M = "Reads/unique UMI",
  
  Preseq_species_richness_estimation_All = "Species richness (all)",
  Preseq_species_richness_estimation_lower_CI_All = "Species richness (all, lower CI)",
  Preseq_species_richness_estimation_upper_CI_All = "Species richness (all, upper CI)",
  Preseq_species_richness_estimation_Tissue = "Species richness (tissue)",
  Preseq_species_richness_estimation_lower_CI_Tissue = "Species richness (tissue, lower CI)",
  Preseq_species_richness_estimation_upper_CI_Tissue = "Species richness (tissue, upper CI)",
  Preseq_species_richness_estimation_per_mm2_tissue = "Species richness/mm² (tissue)",
  Preseq_species_richness_estimation_lower_CI_per_mm2_tissue = "Species richness/mm² (tissue, lower CI)",
  Preseq_species_richness_estimation_upper_CI_per_mm2_tissue = "Species richness/mm² (tissue, upper CI)"
  
)


# Define colors
tech_colors <- c(
  "BMKMANU S1000" = "#000000",
  "BMKMANU S3000" = "#CC79A7",
  "Curio Seeker" = "#9B9B9B",
  "Nova-ST" = "#D55E00",
  "Stereo-Seq" = "#0072B2",
  "Visium (FFPE)" = "#E69F00",
  "Visium HD (FFPE)" = "#009E73",
  "Visium HD (FF)" = "#F0E442",
  "Visium HD 3'" = "#56B4E9"
)

# Function to create bar plots
make_bar_plot <- function(df, ycol) {
  plot_title <- if (ycol %in% names(pretty_names)) pretty_names[[ycol]] else ycol
  
  # Create df_plot
  df_plot <- df %>%
    select(Sample, Technology, Replicate, !!sym(ycol)) %>%
    rename(y_val = !!sym(ycol)) %>%
    mutate(
      Sample = factor(Sample),
      Sample_num = as.numeric(Sample),
      label_text = case_when(
        is.na(y_val) ~ "NA",
        y_val > 1000 ~ scales::comma(y_val),
        TRUE ~ as.character(y_val)
      )
    )
  
    shade_df <- df_plot %>%
      distinct(Sample, Sample_num) %>%
      mutate(
        grey_bg = Sample %in% c("Mouse brain", "Mouse lung", "Rice root")
      )

  # Specify dodge width
  dodge_width <- 0.9
  
  # Determine if Percentage
  percentage = grepl("Percentage|Sequencing_saturation", ycol)
  
  # Determine y-axis formatting
  max_val <- max(df_plot$y_val, na.rm = TRUE)
  y_labels_func <- if (max_val > 1e6) {
    scales::label_number(scale = 1e-6, big.mark = "", suffix = "M")
  } else {
    waiver()  # default formatting
  }

  # Plot
  p <- ggplot(df_plot, aes(x = Sample, y = y_val, fill = Technology)) +
    geom_rect(
      data = shade_df,
      aes(
        xmin = Sample_num - 0.5,
        xmax = Sample_num + 0.5,
        ymin = -Inf,
        ymax = Inf,
        fill = NULL,          # don't use Technology colors here
        alpha = grey_bg
      ),
      inherit.aes = FALSE,
      color = NA
    ) +
    scale_alpha_manual(
      values = c(`TRUE` = 0.12, `FALSE` = 0),
      guide = "none"
    ) +
    geom_bar(aes(group = Technology), stat = "identity", position = position_dodge2(width = dodge_width, preserve = "single", padding=0), color = "black", width=0.9) +
    geom_text(aes(label = label_text, group = Technology), position = position_dodge2(width = dodge_width, preserve = "single"), hjust = -0.1, vjust = 0.5, size = 2.6) +
    geom_text(aes(label = Replicate, group = Technology), position = position_dodge2(width = dodge_width, preserve = "single"), y = 0.005 * max(df_plot$y_val, na.rm = TRUE), hjust = 0, vjust = 0.5, size = 2.6) +
    coord_flip(clip = "off") +
    theme_minimal() +
    theme(
      axis.line.x = element_line(linewidth = 0.4),
      axis.line.y = element_line(linewidth = 0.4),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(color = "black", size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey30", linewidth = 0.25),
      axis.text = element_text(color = "black", size = 11),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    labs(title = NULL, y = str_wrap(plot_title, width = 35), x = NULL, fill = "Technology") +
    scale_fill_manual(values = tech_colors) +
    guides(fill = guide_legend(reverse = TRUE))
  
  p <- p + scale_y_continuous(
    breaks = if (percentage) seq(0, 101, by = 10) else pretty_breaks(n = 8),
    labels = y_labels_func,
    limits = if (percentage) c(0, 101) else NULL,
    expand = expansion(mult = if (percentage) c(0, 0.1) else c(0, 0.25))
  )
  
  ggsave(paste0(gsub(" ", "_", ycol), ".png"), plot = p, width = 8, height = 7)
}

# Loop over all columns
for (col_name in cols_to_plot) {
  make_bar_plot(df, col_name)
}
