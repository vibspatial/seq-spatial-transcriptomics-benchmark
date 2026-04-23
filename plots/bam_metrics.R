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
  ~Sample,      ~sample_type,    ~Technology,       ~replicate,    ~bam_dir,
  "SCA013",     "Mouse_lung",     "Visium",         NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/bam_output",
  "SCA014",     "Mouse_liver",    "Visium",         NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/bam_output",
  "SCA016",     "Mouse_brain",    "Visium",         NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA016/bam_output",
  "SCA015",     "Human_lung",     "Visium",         NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/bam_output",
  
  "SCA001",     "Mouse_lung",     "Visium_HD_FFPE", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA001/bam_output",
  "SCA002",     "Mouse_liver",    "Visium_HD_FFPE", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA002/bam_output",
  "SCA003",     "Mouse_brain",    "Visium_HD_FFPE", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA003/bam_output",
  #"SCA005",     "Mouse_colon",    "Visium_HD_FFPE", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA005/bam_output",
  "SCA010",     "Human_lung",     "Visium_HD_FFPE", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010/bam_output",
  
  "SCA017",     "Mouse_lung",     "Visium_HD_FF",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA017/bam_output",
  "SCA018",     "Mouse_liver",    "Visium_HD_FF",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA018/bam_output",
  "SCA019",     "Mouse_brain",    "Visium_HD_FF",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA019/bam_output",
  "SCA020",     "Human_lung",     "Visium_HD_FF",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA020/bam_output",
  
  "SCA021",     "Mouse_liver",    "Visium_HD_3",    1,              "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA021/bam_output",
  "SCA023",     "Mouse_liver",    "Visium_HD_3",    2,              "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA023/bam_output",
  "SCA024",     "Mouse_brain",    "Visium_HD_3",    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA024/bam_output",
  "FPE010",     "Rice_root",      "Visium_HD_3",    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/FPE010/bam_output",
  
  "Exp74",      "Mouse_brain",    "BMKMANU_S1000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/exp74/bam_output",
  
  "Exp101",     "Mouse_lung",     "BMKMANU_S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/bam_output",
  "Exp100",     "Mouse_liver",    "BMKMANU_S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/bam_output",
  "Exp99",      "Mouse_brain",    "BMKMANU_S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/bam_output",
  "Exp102",     "Human_lung",     "BMKMANU_S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/bam_output",
  "psb",        "Rice_root",      "BMKMANU_S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/psb/bam_output",
  
  "SC12",       "Mouse_lung",     "Curio_Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC12/bam_output",
  "SC13",       "Mouse_liver",    "Curio_Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC13/bam_output",
  "SC11",       "Mouse_brain",    "Curio_Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/bam_output",
  "SC14",       "Human_lung",     "Curio_Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC14/bam_output",
  
  "Exp93-sampleSPC002",     "Mouse_lung",     "Nova-ST",        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC002/bam_output",
  "Exp100-sampleSPC004",    "Mouse_liver",    "Nova-ST",        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp100-sampleSPC004/bam_output",
  "Exp65-sampleSPC014",     "Mouse_brain",    "Nova-ST",        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp65-sampleSPC014/bam_output",
  "Exp93-sampleSPC022",     "Human_lung",     "Nova-ST",        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC022/bam_output",
  "Exp100-samplePSB",       "Rice_root",      "Nova-ST",        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp100-samplePSB/bam_output",
  
  "spc002",     "Mouse_lung",     "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc002/bam_output",
  "spc004",     "Mouse_liver",    "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc004/bam_output",
  "spc014",     "Mouse_brain",    "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc014/bam_output",
  "spc022",     "Human_lung",     "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc022/bam_output",
  "psb_rice",   "Rice_root",      "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/psb/bam_output",
  
)

# Columns ------------------------------------------------------
read_columns <- c(
  "number_of_spots_with_at_least_1_UMI_in_tissue" = "Spots in tissue with ≥1 UMI count",
  "percentage_spots_recovered" = "Percentage of spots recovered in tissue (%)",
  "total_number_of_genes" = "Genes",
  "total_number_of_UMIs" = "UMI counts",
  "total_number_of_unique_UMIs" = "Unique UMI counts",
  "number_of_genes_in_tissue" = "Genes detected in tissue",
  "number_of_unique_UMIs_in_tissue" = "Unique UMI counts detected in tissue",
  "percentage_of_genes_in_tissue" = "Percentage of genes detected in tissue (%)",
  "percentage_of_unique_UMIs_in_tissue" = "Percentage of unique UMI counts detected in tissue (%)",
  "sequencing_saturation_reads" = "Sequencing saturation (%)",
  "sequencing_saturation_UMIs" = "Sequencing saturation (%)",
  "unique_UMIs_per_mm2_in_tissue" = "Unique UMI counts detected in tissue/mm²",
  "median_UMIs_per_gene" = "Median UMI counts per gene",
  "median_UMIs_per_gene_per_mm2" = "Median UMI counts per gene/mm²",
  "median_unique_UMIs_per_gene" = "Median unique UMI counts per gene",
  "median_unique_UMIs_per_gene_per_mm2" = "Median unique UMI counts per gene/mm²",
  "median_unique_UMIs_per_gene_in_tissue" = "Median unique UMI counts in tissue per gene",
  "median_unique_UMIs_per_gene_in_tissue_per_mm2" = "Median unique UMI counts in tissue per gene/mm²",
  "percentage_of_reads_with_usable_UMI" = "Percentage of reads contributing to UMI counts (%)",
  "ratio_reads_per_unique_UMI" = "Reads/unique UMI count",
  "median_reads_per_unique_UMI" = "Median reads/unique UMI count",
  "median_reads_per_unique_UMI_tissue" = "Median reads/unique UMI count in tissue"
)

umi_columns <- c(
  "total_number_of_genes" = "Genes",
  "total_number_of_unique_UMIs" = "Unique UMI counts",
  "number_of_unique_UMIs_in_tissue" = "Unique UMI counts detected in tissue",
  "unique_UMIs_per_mm2_in_tissue" = "Unique UMI counts detected in tissue/mm²",
  "sequencing_saturation_UMIs" = "Sequencing saturation (%)",
  "median_unique_UMIs_per_gene_in_tissue_per_mm2" = "Median unique UMI counts in tissue per gene/mm²"
)

# Colors -------------------------------------------------------
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

sample_type_colors <- c(
  "Mouse lung" = "#56B4E9",
  "Mouse liver" = "#D55E00",
  "Mouse brain" = "#CC79A7",
  "Human lung" = "#0072B2",
  "Rice root" = "#009E73"
)


# Output directories -----------------------------------------------------------------
tech_outdirs <- list(
  "Visium (FFPE)" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/plots/bam_plots",
  "Visium HD (FFPE)" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/plots/FFPE/bam_plots",
  "Visium HD (FF)" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/plots/FF/bam_plots",
  "Visium HD 3'" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/plots/3_prime/bam_plots",
  "BMKMANU S1000" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/plots/bam_plots",
  "BMKMANU S3000" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/plots/bam_plots",
  "Curio Seeker" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/plots/bam_plots",
  "Nova-ST" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/plots/bam_plots",
  "Stereo-Seq" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/plots/bam_plots"
)

# Define function to read all metrics.csv for a given sample -------------------------
read_sample_metrics <- function(bam_dir, sample_name, sample_type, technology, replicate) {
  files <- list.files(file.path(bam_dir, "output_metrics"),
                      pattern = "metrics.csv", recursive = TRUE, full.names = TRUE)
  
  df <- map_dfr(files, ~{
    read_csv(.x, show_col_types = FALSE) %>%
      mutate(
        target_reads = str_extract(.x, "sample_size_\\d+") %>%
          str_remove("sample_size_") %>%
          as.numeric()
      )
  })
  
  df$Sample <- sample_name
  df$sample_type <- sample_type
  df$Technology <- technology
  df$bam_dir <- bam_dir
  df$Replicate <- replicate
  return(df)
}

# Load all samples --------------------------------------------------------
all_data <- map_dfr(1:nrow(sample_metadata), function(i) {
  row <- sample_metadata[i, ]
  read_sample_metrics(row$bam_dir, row$Sample, row$sample_type, row$Technology, row$replicate)
})

zero_rows <- all_data %>%
  group_by(Sample, sample_type, Technology, bam_dir) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    across(where(is.numeric) & !matches("Replicate"), ~0),
    total_number_of_reads = 0,
    Replicate = NA
  )

all_data <- bind_rows(all_data, zero_rows) %>%
  mutate(Technology = recode(Technology,
                             "Visium_HD_3" = "Visium HD 3'",
                             "Visium_HD_FF" = "Visium HD (FF)",
                             "Visium_HD_FFPE" = "Visium HD (FFPE)",
                             "Visium" = "Visium (FFPE)",
                             "Stereo-Seq" = "Stereo-Seq",
                             "Nova-ST" = "Nova-ST",
                             "Curio_Seeker" = "Curio Seeker",
                             "BMKMANU_S3000" = "BMKMANU S3000",
                             "BMKMANU_S1000" = "BMKMANU S1000"),
         Technology = factor(Technology, levels = c(
           "Visium HD 3'",
           "Visium HD (FF)",
           "Visium HD (FFPE)",
           "Visium (FFPE)",
           "Stereo-Seq",
           "Nova-ST",
           "Curio Seeker",
           "BMKMANU S3000",
           "BMKMANU S1000")),
         sample_type = recode(sample_type,
                              "Rice_root" = "Rice root",
                              "Human_lung" = "Human lung",
                              "Mouse_lung" = "Mouse lung",
                              "Mouse_liver" = "Mouse liver",
                              "Mouse_brain" = "Mouse brain"),
         sample_type = factor(sample_type, levels = c(
           "Rice root",
           "Human lung",
           "Mouse lung",
           "Mouse liver",
           "Mouse brain")),
         `Sample type` = sample_type)

# Function to make line plots ---------------------------------------------
make_line_plot <- function(df, xcol, ycol, outpath, color_by, palette) {
  df_plot <- df %>%
    select(Sample, `Sample type`, Technology, !!sym(xcol), !!sym(ycol)) %>%
    rename(x_val = !!sym(xcol)) %>%
    rename(y_val = !!sym(ycol))
  
  percentage <- grepl(
    "percentage|sequencing_saturation", ycol)
  max_val <- max(df_plot$y_val, na.rm = TRUE)
  y_labels_func <- if (max_val > 1e6) {
    scales::label_number(scale = 1e-6, big.mark = "", suffix = "M")
  } else {
    scales::comma
  }
  
  y_label <- read_columns[[ycol]]
  if (is.null(y_label)) y_label <- umi_columns[[ycol]]
  if (is.null(y_label)) y_label <- ycol
  
  if (xcol == "total_number_of_reads") x_label <- "Reads" else xcol
  if (xcol == "total_number_of_UMIs") x_label <- "UMI counts" else xcol
  
  # Visium HD 3' liver labels
  label_df <- df_plot %>%
    left_join(df %>% select(Sample, Replicate) %>% distinct(), by = "Sample") %>%
    filter(!is.na(Replicate)) %>%
    group_by(Sample) %>%
    filter(x_val == max(x_val, na.rm = TRUE)) %>%
    ungroup()
  
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
      axis.text = element_text(color = "black", size = 11),
      axis.title.x=element_text(size = 12, margin = margin(r = 5)),
      axis.title.y=element_text(size = 12),
      plot.margin = margin(10, 20, 10, 10),
      #legend.position = "none"
    ) +
    labs(title = NULL, x = x_label, y = str_wrap(y_label, width = 35), color = color_by) +
    scale_color_manual(values = palette) +
    scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, big.mark = "", suffix = "M"),
      breaks = c(0, 200e6, 400e6, 600e6, 800e6, 1000e6),
      limits = c(0, 1000e6),
      expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(
      breaks = if (percentage) seq(0, 100, by = 10) else pretty_breaks(n = 6),
      labels = y_labels_func,
      limits = if (percentage) c(0, 100) else c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    guides(color = guide_legend(reverse = TRUE))
  
    p <- p +
      geom_text(
        data = label_df,
        aes(label = Replicate),
        hjust = -0.2,
        size = 3,
        show.legend = FALSE,
        color = "black", 
      )
  
  dir.create(file.path(outpath, xcol), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outpath, xcol, paste0(gsub(" ", "_", ycol), ".png")),
         plot = p, width = 6, height = 4)
}

# Generate plots per sample -----------------------------------------------
samples <- distinct(all_data, Sample, bam_dir)
for (i in seq_len(nrow(samples))) {
  sub <- filter(all_data, Sample == samples$Sample[i])
  outdir <- file.path(samples$bam_dir[i], "plots")
  for (col in names(read_columns)) {
    make_line_plot(sub, "total_number_of_reads", col, outdir, "Technology", tech_colors)
  }
  for (col in names(umi_columns)) {
    make_line_plot(sub, "total_number_of_UMIs", col, outdir, "Technology", tech_colors)
  }
}

# Generate plots per technology -------------------------------------------
tech_list <- unique(all_data$Technology)
for (tech in tech_list) {
  subset_df <- all_data[all_data$Technology == tech, ]
  
  outdir <- outdir <- tech_outdirs[[tech]]
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (col in names(read_columns)) {
    make_line_plot(subset_df, "total_number_of_reads", col, outdir, "Sample type", sample_type_colors)
  }
  for (col in names(umi_columns)) {
    make_line_plot(subset_df, "total_number_of_UMIs", col, outdir, "Sample type", sample_type_colors)
  }
}

# Generate plots per sample type ------------------------------------------
sample_type_list <- unique(all_data$sample_type)
for (sample_type in sample_type_list) {
  subset_df <- all_data[all_data$sample_type == sample_type, ]
  
  outdir <- file.path("/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/plots/bam_plots", sample_type)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (col in names(read_columns)) {
    make_line_plot(subset_df, "total_number_of_reads", col, outdir, "Technology", tech_colors)
  }
  for (col in names(umi_columns)) {
    make_line_plot(subset_df, "total_number_of_UMIs", col, outdir, "Technology", tech_colors)
  }
}