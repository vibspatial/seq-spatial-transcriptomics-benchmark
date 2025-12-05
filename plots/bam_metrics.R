# Working directory -------------------------------------------------------
setwd("E:/2024-07-UTBenchmark-SpC")

# Import libraries --------------------------------------------------------
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(scales)

# Define look-up table ------------------------------------------
sample_metadata <- tribble(
  ~Sample,      ~sample_type,    ~Technology,       ~bam_dir,                                                                    ~effective_capture_area_in_tissue,
  "SCA013",     "Mouse_lung",     "Visium",         "E:/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/bam_output",        4.2,
  "SCA014",     "Mouse_liver",    "Visium",         "E:/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/bam_output",        7.6,
  "SCA016",     "Mouse_brain",    "Visium",         "E:/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA016/bam_output",        6.6,
  "SCA015",     "Human_lung",     "Visium",         "E:/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/bam_output",        9.3,
  
  "SCA001",     "Mouse_lung",     "Visium_HD_FFPE", "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA001/bam_output",     15.3,
  "SCA002",     "Mouse_liver",    "Visium_HD_FFPE", "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA002/bam_output",     29.4,
  "SCA003",     "Mouse_brain",    "Visium_HD_FFPE", "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA003/bam_output",     24.9,
  "SCA005",     "Mouse_colon",    "Visium_HD_FFPE", "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA005/bam_output",     7,
  "SCA010",     "Human_lung",     "Visium_HD_FFPE", "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010/bam_output",     33.3,
  
  "SCA017",     "Mouse_lung",     "Visium_HD_FF",   "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA017/bam_output",     10.4,   
  "SCA018",     "Mouse_liver",    "Visium_HD_FF",   "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA018/bam_output",     21.2,
  "SCA019",     "Mouse_brain",    "Visium_HD_FF",   "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA019/bam_output",     19.1,
  "SCA020",     "Human_lung",     "Visium_HD_FF",   "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA020/bam_output",     40.1,
  
  "Exp74",      "Mouse_brain",    "BMKMANU_S1000",  "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/Exp74/bam_output",      4.7,
  
  "Exp101",     "Mouse_lung",     "BMKMANU_S3000",  "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/bam_output",     4.7,
  "Exp100",     "Mouse_liver",    "BMKMANU_S3000",  "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/bam_output",     6.8,
  "Exp99",      "Mouse_brain",    "BMKMANU_S3000",  "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/bam_output",      7.9,
  "Exp102",     "Human_lung",     "BMKMANU_S3000",  "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/bam_output",     9.6,
  
  "SC12",       "Mouse_lung",     "Curio_Seeker",   "E:/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC12/bam_output",    5.5,
  "SC13",       "Mouse_liver",    "Curio_Seeker",   "E:/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC13/bam_output",    5.4,
  "SC11",       "Mouse_brain",    "Curio_Seeker",   "E:/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/bam_output",    5.5,
  "SC14",       "Human_lung",     "Curio_Seeker",   "E:/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC14/bam_output",    5.6,
)

# Columns ------------------------------------------------------
read_columns <- c(
  "number_of_spots_with_at_least_1_UMI_in_tissue",
  "percentage_spots_recovered",
  "total_number_of_genes",
  "total_number_of_UMIs",
  "total_number_of_unique_UMIs",
  "number_of_genes_in_tissue",
  "number_of_unique_UMIs_in_tissue",
  "percentage_of_genes_in_tissue",
  "percentage_of_unique_UMIs_in_tissue",
  "sequencing_saturation_reads",
  "sequencing_saturation_UMIs",
  "unique_UMIs_per_mm2_in_tissue",
  "median_UMIs_per_gene",
  "median_UMIs_per_gene_per_mm2",
  "median_unique_UMIs_per_gene",
  "median_unique_UMIs_per_gene_per_mm2",
  "median_unique_UMIs_per_gene_in_tissue",
  "median_unique_UMIs_per_gene_in_tissue_per_mm2",
  "percentage_of_reads_with_usable_UMI",
  "ratio_reads_per_unique_UMI",
  "median_reads_per_unique_UMI",
  "median_reads_per_unique_UMI_tissue",
  "unique_UMIs_per_mm2_in_tissue_corrected"
)

umi_columns <- c(
  "total_number_of_genes",
  "total_number_of_unique_UMIs",
  "number_of_unique_UMIs_in_tissue",
  "unique_UMIs_per_mm2_in_tissue",
  "sequencing_saturation_UMIs",
  "median_unique_UMIs_per_gene_in_tissue_per_mm2",
  "unique_UMIs_per_mm2_in_tissue_corrected"
)

# Colors -------------------------------------------------------
tech_colors <- c(
  "Visium" = "#F4B18C",
  "Visium_HD_FFPE" = "#A6E18F",
  "Visium_HD_FF" = "#FFEC61",
  "BMKMANU_S1000" = "#94BEE8",
  "BMKMANU_S3000" = "#DC80D3",
  "Curio_Seeker" = "#9B9B9B"
)

# consistent sample_type colors across technologies
sample_type_colors <- c(
  "Mouse_lung" = "#999999",
  "Mouse_liver" = "#E69F00",
  "Mouse_brain" = "#56B4E9",
  "Human_lung" = "#009E73",
  "Mouse_colon" = "#F0E442" # "#0072B2", "#D55E00", "#CC79A7"
)

# Output directories -----------------------------------------------------------------
tech_outdirs <- list(
  "Visium" = "E:/2024-07-UTBenchmark-SpC/data/processed/Visium/plots/bam_plots",
  "Visium_HD_FFPE" = "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/plots/FFPE/bam_plots",
  "Visium_HD_FF" = "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/plots/FF/bam_plots",
  "BMKMANU_S1000" = "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/plots/bam_plots",
  "BMKMANU_S3000" = "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/plots/bam_plots",
  "Curio_Seeker" = "E:/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/plots/bam_plots"
)

# Define function to read all metrics.csv for a given sample -------------------------
read_sample_metrics <- function(bam_dir, sample_name, sample_type, technology, effective_capture_area_in_tissue) {
  files <- list.files(file.path(bam_dir, "output_metrics"),
                      pattern = "metrics.csv", recursive = TRUE, full.names = TRUE)
  
  df <- map_dfr(files, ~{
    read_csv(.x, show_col_types = FALSE) %>%
      mutate(
        target_reads = str_extract(.x, "sample_size_\\d+") %>%
          str_remove("sample_size_") %>%
          as.numeric()
      ) %>%
      mutate(
        unique_UMIs_per_mm2_in_tissue_corrected = !!sym("number_of_unique_UMIs_in_tissue") / effective_capture_area_in_tissue
      )
  })
  
  df$Sample <- sample_name
  df$sample_type <- sample_type
  df$Technology <- technology
  df$bam_dir <- bam_dir
  return(df)
}

# Load all samples --------------------------------------------------------
all_data <- map_dfr(1:nrow(sample_metadata), function(i) {
  row <- sample_metadata[i, ]
  read_sample_metrics(row$bam_dir, row$Sample, row$sample_type, row$Technology, row$effective_capture_area_in_tissue)
})

zero_rows <- all_data %>%
  group_by(Sample, sample_type, Technology, bam_dir) %>%
  slice(1) %>%                          # copy structure
  ungroup() %>%
  mutate(across(where(is.numeric), ~0), # set all numeric cols to 0
         total_number_of_reads = 0)

all_data <- bind_rows(all_data, zero_rows)

# Function to make line plots ---------------------------------------------
make_line_plot <- function(df, xcol, ycol, outpath, color_by, palette) {
  df_plot <- df %>%
    select(Sample, sample_type, Technology, !!sym(xcol), !!sym(ycol)) %>%
    rename(x_val = !!sym(xcol)) %>%
    rename(y_val = !!sym(ycol))
  
  percentage <- grepl(
    "percentage|sequencing_saturation", ycol)
  max_val <- max(df_plot$y_val, na.rm = TRUE)
  y_labels_func <- if (max_val > 1e6) {
    scales::label_number(scale = 1e-6, suffix = "M")
  } else {
    scales::comma
  }
  
  p <- ggplot(df_plot, aes(x = x_val, y = y_val,
                           color = !!sym(color_by), group = Sample)) +
    geom_line(size=1.1) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      text = element_text(color = "black", size = 12),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(size=0.2, color="#DFDFDF"),
      panel.grid.major.y = element_line(size=0.2, color="#DFDFDF"),
      axis.text = element_text(color = "black", size = 12),
      axis.title.x=element_text(size = 12, margin = margin(r = 5)),
      axis.title.y=element_blank(),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    labs(title = ycol, x = xcol, y = ycol, color = color_by) +
    scale_color_manual(values = palette) +
    scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = "M"),
      breaks = pretty_breaks(n = 6),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(
      breaks = if (percentage) seq(0, 100, by = 10) else pretty_breaks(n = 6),
      labels = y_labels_func,
      limits = if (percentage) c(0, 100) else c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    guides(color = guide_legend(reverse = TRUE))
  
  dir.create(file.path(outpath, xcol), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outpath, xcol, paste0(gsub(" ", "_", ycol), ".png")),
         plot = p, width = 6, height = 4)
}

# Generate plots per sample -----------------------------------------------
samples <- distinct(all_data, Sample, bam_dir)
for (i in seq_len(nrow(samples))) {
  sub <- filter(all_data, Sample == samples$Sample[i])
  outdir <- file.path(samples$bam_dir[i], "plots")
  for (col in read_columns) {
    make_line_plot(sub, "total_number_of_reads", col, outdir, "Technology", tech_colors)
  }
  for (col in umi_columns) {
    make_line_plot(sub, "total_number_of_UMIs", col, outdir, "Technology", tech_colors)
  }
}

# Generate plots per technology -------------------------------------------
tech_list <- unique(all_data$Technology)
for (tech in tech_list) {
  subset_df <- all_data[all_data$Technology == tech, ]
  
  outdir <- outdir <- tech_outdirs[[tech]]
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (col in read_columns) {
    make_line_plot(subset_df, "total_number_of_reads", col, outdir, "sample_type", sample_type_colors)
  }
  for (col in umi_columns) {
    make_line_plot(subset_df, "total_number_of_UMIs", col, outdir, "sample_type", sample_type_colors)
  }
}

# Generate plots per sample type ------------------------------------------
sample_type_list <- unique(all_data$sample_type)
for (sample_type in sample_type_list) {
  subset_df <- all_data[all_data$sample_type == sample_type, ]
  
  outdir <- file.path("E:/2024-07-UTBenchmark-SpC/data/processed/plots/bam_plots", sample_type)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (col in read_columns) {
    make_line_plot(subset_df, "total_number_of_reads", col, outdir, "Technology", tech_colors)
  }
  for (col in umi_columns) {
    make_line_plot(subset_df, "total_number_of_UMIs", col, outdir, "Technology", tech_colors)
  }
}

# Bar plots ---------------------------------------------------
make_bar_plot <- function(df, ycol, outpath) {
  df_plot <- df %>%
    group_by(Sample) %>%
    filter(total_number_of_reads == max(total_number_of_reads, na.rm = TRUE)) %>%
    ungroup() %>%
    select(Sample, Technology, sample_type, !!sym(ycol)) %>%
    rename(y_val = !!sym(ycol)) %>%
    mutate(
      label_text = case_when(
        is.na(y_val) ~ "NA",
        TRUE ~ sprintf("%.1f", y_val)
      )
    )
  
  # Make sample_type a factor with the same order as the color vector
  df_plot <- df_plot %>%
    mutate(sample_type = factor(sample_type, levels = rev(names(sample_type_colors))))
  
  # Reorder Sample by sample_type (so bars follow the legend order)
  df_plot <- df_plot %>%
    mutate(Sample = forcats::fct_reorder(Sample, as.numeric(sample_type)))
  
  # Specify dodge width
  dodge_width <- 0.9
  
  # Plot
  p <- ggplot(df_plot, aes(x = Sample, y = y_val, fill = sample_type)) +
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
    labs(title = ycol, y = NULL, x = NULL, fill = "sample_type") +
    scale_fill_manual(values = sample_type_colors) +
    guides(fill = guide_legend(reverse = TRUE))
  
  p <- p + scale_y_continuous(
    breaks = seq(0, 101, by = 10),
    limits = c(0, 101),
    expand = expansion(mult = c(0, 0.05))
  )
  
  dir.create(file.path(outpath, "bar_plots"), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outpath, "bar_plots", paste0(ycol, ".png")),
         plot = p, width = 8, height = 2.5)
}

# Curio
curio_cols <- c("reads_Assigned",
                "reads_Unassigned_NoFeatures",
                "reads_Unassigned_MultiMapping",
                "reads_Unassigned_Ambiguity",
                "reads_NA", 
                "percentage_of_genes_in_tissue",
                "percentage_of_reads_with_usable_UMI",
                "percentage_of_unique_UMIs_in_tissue",
                "percentage_spots_recovered")

curio_df <- all_data %>% filter(Technology == "Curio_Seeker")
curio_outdir <- tech_outdirs[["Curio_Seeker"]]
for (col_name in curio_cols) {
  make_bar_plot(curio_df, col_name, curio_outdir)
}

# Visium
visium_cols <- c("reads_mapped",
                 "reads_mapped_confident",
                 "reads_mapped_confident_to_filtered",
                 "reads_half_mapped",
                 "reads_split_mapped",
                 "Q30_bases_read",
                 "Q30_bases_spatial_barcode",
                 "Q30_bases_UMI",
                 "percentage_of_genes_in_tissue",
                 "percentage_of_reads_with_usable_UMI",
                 "percentage_of_unique_UMIs_in_tissue",
                 "percentage_spots_recovered")

for (tech in c("Visium", "Visium_HD_FFPE", "Visium_HD_FF")) {
  subset_df <- all_data %>% filter(Technology == tech)
  outdir <- tech_outdirs[[tech]]
  for (col_name in visium_cols) {
    make_bar_plot(subset_df, col_name, outdir)
  }
}