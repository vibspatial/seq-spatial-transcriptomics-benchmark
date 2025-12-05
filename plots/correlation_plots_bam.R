# Libraries ---------------------------------------------------------------
library(lattice)
library(latticeExtra)
library(tidyverse)
library(stringi)
library(data.table)
library(ggplot2)
library(openxlsx)
library(GGally)
library(viridis)
library(ggpubr)

# Root directories --------------------------------------------------------
root_dir <- "E:/2024-07-UTBenchmark-SpC/data/processed"
out_root <- file.path(root_dir, "plots/bam_plots")

# Define look-up table ------------------------------------------
sample_metadata <- tribble(
  ~Sample,      ~sample_type,    ~Technology,       ~bam_dir,
  "SCA013",     "Mouse_lung",     "Visium",         "E:/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/bam_output",
  #"SCA014",     "Mouse_liver",    "Visium",         "E:/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/bam_output",
  "SCA016",     "Mouse_brain",    "Visium",         "E:/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA016/bam_output",
  "SCA015",     "Human_lung",     "Visium",         "E:/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/bam_output",
  
  "SCA001",     "Mouse_lung",     "Visium_HD_FFPE", "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA001/bam_output",
  "SCA002",     "Mouse_liver",    "Visium_HD_FFPE", "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA002/bam_output",
  "SCA003",     "Mouse_brain",    "Visium_HD_FFPE", "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA003/bam_output",
  "SCA010",     "Human_lung",     "Visium_HD_FFPE", "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010/bam_output",
  
  "SCA017",     "Mouse_lung",     "Visium_HD_FF",   "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA017/bam_output",
  "SCA018",     "Mouse_liver",    "Visium_HD_FF",   "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA018/bam_output",
  "SCA019",     "Mouse_brain",    "Visium_HD_FF",   "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA019/bam_output",
  "SCA020",     "Human_lung",     "Visium_HD_FF",   "E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA020/bam_output",
  
  "Exp74",      "Mouse_brain",    "BMKMANU_S1000",  "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/Exp74/bam_output",
  
  "Exp101",     "Mouse_lung",     "BMKMANU_S3000",  "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/bam_output",
  #"Exp100",     "Mouse_liver",    "BMKMANU_S3000",  "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/bam_output",
  "Exp99",      "Mouse_brain",    "BMKMANU_S3000",  "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/bam_output",
  "Exp102",     "Human_lung",     "BMKMANU_S3000",  "E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/bam_output",
  
  "SC12",       "Mouse_lung",     "Curio_Seeker",   "E:/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC12/bam_output",
  "SC13",       "Mouse_liver",    "Curio_Seeker",   "E:/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC13/bam_output",
  "SC11",       "Mouse_brain",    "Curio_Seeker",   "E:/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/bam_output",
  "SC14",       "Human_lung",     "Curio_Seeker",   "E:/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC14/bam_output",
)

# Read gene_counts per technology and subsample ---------------
read_gene_counts <- function(technology, bam_dir, subsample) {
  path <- file.path(bam_dir, "output_metrics", subsample, "gene_counts.csv")
  if (file.exists(path)) {
    df <- read.csv(path) %>%
      select(gene, unique_gene_counts_tissue_corrected) %>%
      rename(!!technology := unique_gene_counts_tissue_corrected)
    return(df)
  } else {
    return(NULL)
  }
}

prepare_correlation_data <- function(sample_type, subsample, technologies, sample_metadata_subset) {
  dfs <- list()
  
  # Load all available files
  for (technology in technologies) {
    bam_dir <- sample_metadata_subset$bam_dir[sample_metadata_subset$Technology == technology][1]
    df <- read_gene_counts(technology, bam_dir, subsample)
    
    if (is.null(df)) {
      # Create placeholder if file doesn't exist
      df <- tibble(gene = character(0), !!technology := numeric(0))
    }
    
    df <- df %>%
      add_row(gene = "dummy_gene_1", !!technology := 1) %>%
      add_row(gene = "dummy_gene_2", !!technology := 2)
    
    dfs[[technology]] <- df
  }
  
  # Merge all data frames by 'gene'
  merged <- reduce(dfs, full_join, by = "gene")
  
  # Convert columns to numeric and replace zeros or negatives with NA
  merged <- merged %>%
    mutate(across(all_of(technologies), as.numeric)) %>%
    mutate(across(all_of(technologies), ~ifelse(. <= 0, NA, .)))
  
  return(merged)
}

# Log correlation function ---------------------------------------------
diagfun <- function(data, mapping) {
  colname <- rlang::as_label(mapping$x)
  
  if (sum(!is.na(data[[colname]])) <= 2) {
    opacity <- 0
    line_color <- NA
    } else {
    opacity <- 0.5
    line_color <- "black"
    data <- data[!data$gene %in% c("dummy_gene_1", "dummy_gene_2"), ]
    }
  
  ggplot(data = data, mapping = mapping) + 
    geom_density(fill = "#5A2A81", alpha = opacity, color = line_color) + 
    scale_x_log10(limits = c(1, 100000)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          text = element_text(color = "black", size = 6),
          axis.text = element_text(size = 7, color = "black"), 
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
}

lowerfun <- function(data, mapping) {
  colname_x <- rlang::as_label(mapping$x)
  colname_y <- rlang::as_label(mapping$y)
  
  if (sum(!is.na(data[[colname_x]])) <= 2 | sum(!is.na(data[[colname_y]])) <= 2) {
    opacity <- 0
    } else {
    opacity <- 0.6
    data <- data[!data$gene %in% c("dummy_gene_1", "dummy_gene_2"), ]
    }

  ggplot(data = data, mapping = mapping) +
    geom_point(color = "#AC94C0", size = 0.3, alpha = opacity, na.rm = TRUE) +
    scale_x_log10(limits = c(1, 100000)) +
    scale_y_log10(limits = c(1, 100000)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth=0.6) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(color = "black", size = 6),
          axis.text = element_text(size = 7, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
}

upperfun <- function(data, mapping) {
  colname_x <- rlang::as_label(mapping$x)
  colname_y <- rlang::as_label(mapping$y)
  
  if (sum(!is.na(data[[colname_x]])) <= 2 | sum(!is.na(data[[colname_y]])) <= 2) {
    opacity <- 0
  } else {
    opacity <- 0.6
    data <- data[!data$gene %in% c("dummy_gene_1", "dummy_gene_2"), ]
  }
  
  ggplot(data = data, mapping = mapping) +
    stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", contour = TRUE, alpha = opacity, na.rm = TRUE, n = 500, bins = 30) +
    scale_fill_viridis_c() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth=0.6) +
    scale_x_log10(limits = c(1, 100000)) +
    scale_y_log10(limits = c(1, 100000)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(color = "black", size = 6),
          axis.text = element_text(size = 7, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
    stat_cor(
      aes(label = ..r.label..), 
      method = "pearson",
      label.x.npc = "left", label.y.npc = "top", size = 3, hjust = 0
    )
}

# Main loop ---------------------------------------------------------------
for (sample_type in unique(sample_metadata$sample_type)) {
  message("Processing ", sample_type, " ...")
  
  sample_metadata_subset <- sample_metadata %>%
    filter(sample_type == !!sample_type)
  
  technologies <- unique(sample_metadata_subset$Technology)
  subsample_dirs <- paste0("sample_size_", format(seq(10000000, 1000000000, 10000000), scientific=F, trim = TRUE))

  for (subsample in subsample_dirs) {
    df <- prepare_correlation_data(sample_type, subsample, technologies, sample_metadata_subset)
    
    # Skip if only 1 valid column
    if (sum(colSums(df[technologies], na.rm = TRUE) > 3) <= 1) next
    
    message("  Subsample: ", subsample)
    # Generate correlation plot
    plot_correlation <- ggpairs(
      data = df,
      columns = technologies,
      diag = list(continuous = wrap(diagfun)),
      lower = list(continuous = wrap(lowerfun)),
      upper = list(continuous = wrap(upperfun))
    ) + theme(strip.text = element_text(size = 7))
    
    # Save
    outdir <- file.path(out_root, sample_type, "gene_counts")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    outfile <- file.path(outdir, paste0(sample_type, "_", subsample, "_log_correlation.png"))
    ggsave(outfile, plot_correlation, width = 7, height = 7, dpi = 300)
  }
}
