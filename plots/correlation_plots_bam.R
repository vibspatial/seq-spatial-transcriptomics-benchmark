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
library(ggtext)
library(grid)
library(scales)
library(ggpointdensity)

# Root directories --------------------------------------------------------
root_dir <- "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed"
out_root <- file.path(root_dir, "plots/bam_plots")

# Define look-up table ------------------------------------------
sample_metadata <- tribble(
  ~Sample,      ~sample_type,    ~Technology,       ~replicate,    ~bam_dir,
  "SCA013",     "Mouse lung",     "Visium (FFPE)",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/bam_output",
  "SCA014",     "Mouse liver",    "Visium (FFPE)",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/bam_output",
  "SCA016",     "Mouse brain",    "Visium (FFPE)",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA016/bam_output",
  "SCA015",     "Human lung",     "Visium (FFPE)",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/bam_output",
  
  "SCA001",     "Mouse lung",     "Visium HD (FFPE)", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA001/bam_output",
  "SCA002",     "Mouse liver",    "Visium HD (FFPE)", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA002/bam_output",
  "SCA003",     "Mouse brain",    "Visium HD (FFPE)", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA003/bam_output",
  #"SCA005",     "Mouse colon",    "Visium HD (FFPE)", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA005/bam_output",
  "SCA010",     "Human lung",     "Visium HD (FFPE)", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010/bam_output",
  
  "SCA017",     "Mouse lung",     "Visium HD (FF)",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA017/bam_output",
  "SCA018",     "Mouse liver",    "Visium HD (FF)",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA018/bam_output",
  "SCA019",     "Mouse brain",    "Visium HD (FF)",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA019/bam_output",
  "SCA020",     "Human lung",     "Visium HD (FF)",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA020/bam_output",
  
  "SCA021",     "Mouse liver",    "Visium HD 3'",    1,              "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA021/bam_output",
  "SCA023",     "Mouse liver",    "Visium HD 3'",    2,              "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA023/bam_output",
  "SCA024",     "Mouse brain",    "Visium HD 3'",    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA024/bam_output",
  "FPE010",     "Rice root",      "Visium HD 3'",    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/FPE010/bam_output",
  
  "Exp74",      "Mouse brain",    "BMKMANU S1000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/exp74/bam_output",
  
  "Exp101",     "Mouse lung",     "BMKMANU S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/bam_output",
  "Exp100",     "Mouse liver",    "BMKMANU S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/bam_output",
  "Exp99",      "Mouse brain",    "BMKMANU S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/bam_output",
  "Exp102",     "Human lung",     "BMKMANU S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/bam_output",
  "psb",        "Rice root",      "BMKMANU S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/psb/bam_output",

  "SC12",       "Mouse lung",     "Curio Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC12/bam_output",
  "SC13",       "Mouse liver",    "Curio Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC13/bam_output",
  "SC11",       "Mouse brain",    "Curio Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/bam_output",
  "SC14",       "Human lung",     "Curio Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC14/bam_output",
  
  "Exp93-sampleSPC002",     "Mouse lung",     "Nova-ST",        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC002/bam_output",
  "Exp100-sampleSPC004",    "Mouse liver",    "Nova-ST",        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp100-sampleSPC004/bam_output",
  "Exp65-sampleSPC014",     "Mouse brain",    "Nova-ST",        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp65-sampleSPC014/bam_output",
  "Exp93-sampleSPC022",     "Human lung",     "Nova-ST",        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC022/bam_output",
  "Exp100-samplePSB",       "Rice root",      "Nova-ST",        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp100-samplePSB/bam_output",
  
  "spc002",     "Mouse lung",     "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc002/bam_output",
  "spc004",     "Mouse liver",    "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc004/bam_output",
  "spc014",     "Mouse brain",    "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc014/bam_output",
  "spc022",     "Human lung",     "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc022/bam_output",
  "psb_rice",   "Rice root",      "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/psb/bam_output",
  
)

sample_metadata <- sample_metadata %>%
  mutate(
    Technology_plot = case_when(
      !is.na(replicate) ~ paste0(Technology, " (", replicate, ")"),
      is.na(replicate)  ~ Technology
    )
  )

# Technology colors -------------------------------------------------------
tech_colors <- c(
  "BMKMANU S1000"    = "#000000",
  "BMKMANU S3000"    = "#CC79A7",
  "Curio Seeker"     = "#9B9B9B",
  "Nova-ST"          = "#D55E00",
  "Stereo-Seq"       = "#0072B2",
  "Visium (FFPE)"    = "#E69F00",
  "Visium HD (FFPE)" = "#009E73",
  "Visium HD (FF)"   = "#F0E442",
  "Visium HD 3'"     = "#56B4E9"
)

tech_levels <- c(
  "Visium HD 3'",
  "Visium HD (FF)",
  "Visium HD (FFPE)",
  "Visium (FFPE)",
  "Stereo-Seq",
  "Nova-ST",
  "Curio Seeker",
  "BMKMANU S3000",
  "BMKMANU S1000"
)

# Common genes
common_genes_path <- c(
  "Mouse brain" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/plots/genes/Mouse brain/Mouse brain_common_gene_expression_table.csv",
  "Mouse liver" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/plots/genes/Mouse liver/Mouse liver_common_gene_expression_table.csv",
  "Mouse lung" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/plots/genes/Mouse lung/Mouse lung_common_gene_expression_table.csv",
  "Human lung" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/plots/genes/Human lung/Human lung_common_gene_expression_table.csv",
  "Rice root" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/plots/genes/Rice root/Rice root_common_gene_expression_table.csv"
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
  
  # Read in common genes
  common_genes <- read.csv(common_genes_path[[sample_type]])
  
  # Load all available files
  for (technology in technologies) {
    bam_dir <- sample_metadata_subset$bam_dir[sample_metadata_subset$Technology_plot == technology][1]
    df <- read_gene_counts(technology, bam_dir, subsample)
    
    if (is.null(df)) {
      df <- tibble(gene = character(0), !!technology := numeric(0))
    } else {
      df <- df %>%
        filter(gene %in% common_genes$gene)
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
    mutate(across(all_of(technologies), ~ ifelse(. <= 0, NA, .)))
  
  return(merged)
}

# Calculate fold change
make_pairwise_logfc_table <- function(df, technologies, subsample) {
  
  pair_list <- combn(technologies, 2, simplify = FALSE)
  out <- vector("list", length(pair_list))
  
  for (i in seq_along(pair_list)) {
    tech1 <- pair_list[[i]][1]
    tech2 <- pair_list[[i]][2]
    
    pair_df <- df %>%
      select(gene, all_of(c(tech1, tech2))) %>%
      rename(
        expr_tech1 = all_of(tech1),
        expr_tech2 = all_of(tech2)
      ) %>%
      filter(
        !gene %in% c("dummy_gene_1", "dummy_gene_2"),
        !is.na(expr_tech1),
        !is.na(expr_tech2),
        expr_tech1 > 0,
        expr_tech2 > 0
      ) %>%
      mutate(
        subsample = subsample,
        tech1 = tech1,
        tech2 = tech2,
        FC = expr_tech1 / expr_tech2,
        log2FC = log2(expr_tech1 / expr_tech2),
        avg_expression = (expr_tech1 + expr_tech2) / 2
      ) %>%
      select(
        subsample,
        tech1,
        tech2,
        gene,
        expr_tech1,
        expr_tech2,
        FC,
        log2FC,
        avg_expression
      )
    
    out[[i]] <- pair_df
  }
  
  bind_rows(out)
}
summarise_pairwise_logfc <- function(pairwise_logfc_df) {
  pairwise_logfc_df %>%
    group_by(tech1, tech2) %>%
    summarise(
      median_log2FC = median(log2FC, na.rm = TRUE),
      IQR_log2FC = IQR(log2FC, na.rm = TRUE),
      .groups = "drop"
    )
}

current_pairwise_stats <- NULL

# Plotting correlation function ---------------------------------------------
log_label_fun <- function(x) {
  format(x, scientific = FALSE, trim = TRUE, big.mark = ",")
}

wrap_tech_label <- function(x, width = 10) {
  stringr::str_wrap(x, width = width)
}

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
    scale_x_log10(
      limits = c(1, 100000),
      breaks = c(1, 10, 100, 1000, 10000, 100000),
      labels = log_label_fun,
      minor_breaks = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.15),
      text = element_text(color = "black", size = 6),
      axis.text = element_text(size = 7, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.08, "cm"),
    )
}

blank_panel <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) +
    scale_x_log10(
      limits = c(1, 100000),
      breaks = 10^(0:5),
      labels = log_label_fun,
      minor_breaks = NULL
    ) +
    scale_y_log10(
      limits = c(1, 100000),
      breaks = 10^(0:5),
      labels = log_label_fun,
      minor_breaks = NULL
    ) +
    #geom_abline(
    #  slope = 1, intercept = 0,
    #  linetype = "longdash", color = "black", linewidth = 0.6
    #) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.15),
      panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.15),
      text = element_text(color = "black", size = 6),
      axis.text = element_text(size = 7, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.08, "cm"),
      legend.position = "none")
}

lowerfun <- function(data, mapping) {
  colname_x <- rlang::as_label(mapping$x)
  colname_y <- rlang::as_label(mapping$y)
  
  data <- data[!data$gene %in% c("dummy_gene_1", "dummy_gene_2"), ]
  
  complete_idx <- !is.na(data[[colname_x]]) & !is.na(data[[colname_y]]) &
    data[[colname_x]] > 0 & data[[colname_y]] > 0
  
  if (sum(complete_idx) < 2) {
    return(blank_panel(data, mapping))
  }
  
  data <- data[complete_idx, ]
  
  cor_val <- suppressWarnings(cor(
    data[[colname_x]], data[[colname_y]],
    method = "pearson", use = "complete.obs"
  ))
  cor_label <- sprintf("r = %.2f", cor_val)
  
  ggplot(data = data, mapping = mapping) +
    ggpointdensity::geom_pointdensity(
      size = 0.2,
      adjust = 0.5,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "longdash", color = "black", linewidth = 0.6
    ) +
    annotate(
      "text",
      x = 1.2, y = 90000,
      label = cor_label,
      hjust = 0, vjust = 1,
      size = 2
    ) +
    scale_x_log10(
      limits = c(1, 100000),
      breaks = 10^(0:5),
      labels = log_label_fun,
      minor_breaks = NULL
    ) +
    scale_y_log10(
      limits = c(1, 100000),
      breaks = 10^(0:5),
      labels = log_label_fun,
      minor_breaks = NULL
    ) +
    scale_color_viridis_c(
      option = "D",
      trans = "sqrt",
      name = "Local density"
    ) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.15),
      panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.15),
      text = element_text(color = "black", size = 6),
      axis.text = element_text(size = 7, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.08, "cm"),
      legend.position = "none"
    )
}

upperfun <- function(data, mapping) {
  colname_x <- rlang::as_label(mapping$x)
  colname_y <- rlang::as_label(mapping$y)
  
  if (sum(!is.na(data[[colname_x]])) <= 2 || sum(!is.na(data[[colname_y]])) <= 2) {
    return(blank_panel(data, mapping))
  }
  
  data <- data[!data$gene %in% c("dummy_gene_1", "dummy_gene_2"), ]
  
  complete_idx <- !is.na(data[[colname_x]]) & !is.na(data[[colname_y]]) &
    data[[colname_x]] > 0 & data[[colname_y]] > 0
  
  if (sum(complete_idx) < 2) {
    return(blank_panel(data, mapping))
  }
  
  data <- data[complete_idx, ]
  
  stats_row <- current_pairwise_stats %>%
    filter(tech1 == colname_x, tech2 == colname_y)
  
  if (nrow(stats_row) > 0) {
    log2fc_val <- stats_row$median_log2FC[1]
    iqr_val <- stats_row$IQR_log2FC[1]
  } else {
    stats_row <- current_pairwise_stats %>%
      filter(tech1 == colname_y, tech2 == colname_x)
    
    if (nrow(stats_row) > 0) {
      log2fc_val <- -stats_row$median_log2FC[1]
      iqr_val <- stats_row$IQR_log2FC[1]
    } else {
      log2fc_val <- NA_real_
      iqr_val <- NA_real_
    }
  }
  
  stat_label <- if (!is.na(log2fc_val)) {
    paste0(
      "log2FC=", sprintf("%.1f", log2fc_val),
      "\n",
      "IQR=", sprintf("%.1f", iqr_val)
    )
  } else {
    ""
  }
  
  ggplot(data = data, mapping = mapping) +
    geom_density_2d(
      aes(color = after_stat(level)),
      linewidth = 0.30,
      bins = 15,
      na.rm = TRUE,
      show.legend = TRUE
    ) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "longdash",
      color = "black", linewidth = 0.5
    ) +
    annotate(
      "text",
      x = 1.3, y = 90000,
      label = stat_label,
      hjust = 0, vjust = 1,
      size = 2,
      lineheight = 0.95
    ) +
    scale_x_log10(
      limits = c(1, 100000),
      breaks = 10^(0:5),
      labels = log_label_fun,
      minor_breaks = NULL
    ) +
    scale_y_log10(
      limits = c(1, 100000),
      breaks = 10^(0:5),
      labels = log_label_fun,
      minor_breaks = NULL
    ) +
    scale_color_viridis_c(name = "Density") +
    theme_bw() +
    theme(
      panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.15),
      panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.15),
      text = element_text(color = "black", size = 6),
      axis.text = element_text(size = 7, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.08, "cm"),
      legend.position = "none"
    )
}

# Heatmap helpers ---------------------------------------------------------
make_pairwise_correlation_table <- function(df, technologies, subsample) {
  pair_list <- combn(technologies, 2, simplify = FALSE)
  out <- vector("list", length(pair_list))
  
  for (i in seq_along(pair_list)) {
    tech1 <- pair_list[[i]][1]
    tech2 <- pair_list[[i]][2]
    
    tmp <- df %>%
      select(gene, all_of(c(tech1, tech2))) %>%
      rename(
        expr_tech1 = all_of(tech1),
        expr_tech2 = all_of(tech2)
      ) %>%
      filter(
        !gene %in% c("dummy_gene_1", "dummy_gene_2"),
        !is.na(expr_tech1),
        !is.na(expr_tech2),
        expr_tech1 > 0,
        expr_tech2 > 0
      )
    
    cor_val <- if (nrow(tmp) >= 2) {
      suppressWarnings(cor(
        tmp$expr_tech1,
        tmp$expr_tech2,
        method = "pearson",
        use = "complete.obs"
      ))
    } else {
      NA_real_
    }
    
    out[[i]] <- tibble(
      subsample = subsample,
      tech1 = tech1,
      tech2 = tech2,
      pearson_r = cor_val,
      n_genes_used = nrow(tmp)
    )
  }
  
  bind_rows(out)
}

plot_pairwise_metric_heatmap <- function(pairwise_df, sample_type, subsample, outdir,
                                         metric = c("pearson_r", "median_log2FC"),
                                         tech_order = NULL,
                                         suffix = NULL,
                                         legend_title = NULL) {
  library(tidyverse)
  
  metric <- match.arg(metric)
  
  if (is.null(tech_order)) {
    tech_order <- sort(unique(c(pairwise_df$tech1, pairwise_df$tech2)))
  }
  
  if (is.null(suffix)) {
    suffix <- metric
  }
  
  if (is.null(legend_title)) {
    legend_title <- dplyr::case_when(
      metric == "pearson_r" ~ "Pearson r",
      metric == "median_log2FC" ~ "Median log2FC",
      TRUE ~ metric
    )
  }
  
  base_df <- expand.grid(
    tech1 = tech_order,
    tech2 = tech_order,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()
  
  legend_title_wrapped <- stringr::str_wrap(legend_title, width = 10)
  
  if (metric == "pearson_r") {
    heatmap_df <- base_df %>%
      left_join(
        pairwise_df %>% select(tech1, tech2, pearson_r),
        by = c("tech1", "tech2")
      ) %>%
      left_join(
        pairwise_df %>%
          select(tech1, tech2, pearson_r) %>%
          rename(
            tech1 = tech2,
            tech2 = tech1,
            value_rev = pearson_r
          ),
        by = c("tech1", "tech2")
      ) %>%
      mutate(
        value = coalesce(pearson_r, value_rev),
        value = ifelse(tech1 == tech2, 1, value)
      ) %>%
      select(tech1, tech2, value)
    
    fill_scale <- scale_fill_viridis_c(
      option = "viridis",
      limits = c(0, 1),
      name = legend_title_wrapped
    )
  }
  
  if (metric == "median_log2FC") {
    heatmap_df <- base_df %>%
      left_join(
        pairwise_df %>% select(tech1, tech2, median_log2FC),
        by = c("tech1", "tech2")
      ) %>%
      left_join(
        pairwise_df %>%
          select(tech1, tech2, median_log2FC) %>%
          rename(
            tech1 = tech2,
            tech2 = tech1,
            value_rev = median_log2FC
          ) %>%
          mutate(value_rev = -value_rev),
        by = c("tech1", "tech2")
      ) %>%
      mutate(
        value = coalesce(median_log2FC, value_rev),
        value = ifelse(tech1 == tech2, 0, value)
      ) %>%
      select(tech1, tech2, value)
    
    vmax <- max(abs(heatmap_df$value), na.rm = TRUE)
    
    fill_scale <- scale_fill_gradient2(
      low = "#B40426",
      mid = "white",
      high = "#3B4CC0",
      midpoint = 0,
      limits = c(-6, 6),
      oob = scales::squish,
      name = legend_title_wrapped
    )
  }
  
  heatmap_df <- heatmap_df %>%
    mutate(
      tech1 = factor(tech1, levels = tech_order),
      tech2 = factor(tech2, levels = rev(tech_order))
    )
  
  p <- ggplot(heatmap_df, aes(x = tech1, y = tech2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
    fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(color = "black", size = 12),
      axis.text = element_text(color = "black", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(10, 20, 10, 10)
    )
  
  ggsave(
    file.path(outdir, paste0(sample_type, "_", subsample, "_", suffix, "_heatmap.png")),
    p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

# Main loop ---------------------------------------------------------------
for (sample_type in unique(sample_metadata$sample_type)) {
  message("Processing ", sample_type, " ...")
  
  sample_metadata_subset <- sample_metadata %>%
    filter(sample_type == !!sample_type)
  
  technologies <- sample_metadata_subset %>%
    mutate(
      Technology = factor(Technology, levels = rev(tech_levels))
    ) %>%
    arrange(Technology, Sample) %>%
    pull(Technology_plot)
  
  subsample_dirs <- paste0("sample_size_", format(seq(10000000, 1000000000, 10000000), scientific=F, trim = TRUE))

  for (subsample in subsample_dirs) {
    df <- prepare_correlation_data(sample_type, subsample, technologies, sample_metadata_subset)
    pairwise_logfc_df <- make_pairwise_logfc_table(df, technologies, subsample)
    current_pairwise_stats <- summarise_pairwise_logfc(pairwise_logfc_df)
    pairwise_cor_df <- make_pairwise_correlation_table(df, technologies, subsample)
    
    # Skip if only 1 valid column
    if (sum(colSums(df[technologies], na.rm = TRUE) > 3) <= 1) next
    
    message("  Subsample: ", subsample)
    
    column_labels <- vapply(technologies, wrap_tech_label, character(1), width = 10)
    
    plot_correlation <- ggpairs(
      data = df,
      columns = technologies,
      columnLabels = column_labels,
      diag = list(continuous = wrap(diagfun)),
      lower = list(continuous = wrap(lowerfun)),
      upper = list(continuous = wrap(upperfun))
    ) + 
      theme(
        strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
        strip.text = element_text(size = 7, lineheight = 0.9)
      )
    
    outdir <- file.path(out_root, sample_type, "gene_counts")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    write.csv(
      pairwise_logfc_df,
      file.path(outdir, paste0(sample_type, "_", subsample, "_pairwise_log2FC.csv")),
      row.names = FALSE
    )
    
    write.csv(
      pairwise_cor_df,
      file.path(outdir, paste0(sample_type, "_", subsample, "_pairwise_correlation.csv")),
      row.names = FALSE
    )
    
    outfile <- file.path(outdir, paste0(sample_type, "_", subsample, "_log_correlation.png"))
    ggsave(outfile, plot_correlation, width = 7, height = 7, dpi = 300)
    
    # Heatmaps
    plot_pairwise_metric_heatmap(
      pairwise_df = pairwise_cor_df,
      sample_type = sample_type,
      subsample = subsample,
      outdir = outdir,
      metric = "pearson_r",
      tech_order = technologies,
      suffix = "pearson_r",
      legend_title = "Pearson r"
    )
    
    plot_pairwise_metric_heatmap(
      pairwise_df = current_pairwise_stats,
      sample_type = sample_type,
      subsample = subsample,
      outdir = outdir,
      metric = "median_log2FC",
      tech_order = technologies,
      suffix = "median_log2FC",
      legend_title = "Median log2FC"
    )
  }
}

# Aggregate Log2FC metrics into summary plots
for (sample_type in unique(sample_metadata$sample_type)) {
  message("log2FC over sequencing depth for ", sample_type, " ...")
  
  sample_metadata_subset <- sample_metadata %>%
    filter(sample_type == !!sample_type)
  
  technologies <- sample_metadata_subset %>%
    mutate(
      Technology = factor(Technology, levels = rev(tech_levels))
    ) %>%
    arrange(Technology, Sample) %>%
    pull(Technology_plot)
  
  outdir <- file.path(out_root, sample_type, "gene_counts")
  
  files <- list.files(
    outdir,
    pattern = paste0("^", sample_type, "_sample_size_[0-9]+_pairwise_log2FC\\.csv$"),
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    message("  No pairwise log2FC files found for ", sample_type)
    next
  }
  
  pairwise_logfc_all <- purrr::map_dfr(files, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    
    if (!"log2FC" %in% colnames(df)) {
      return(NULL)
    }
    
    subsample_here <- stringr::str_extract(basename(f), "sample_size_[0-9]+")
    
    df %>%
      mutate(
        sequencing_depth = as.numeric(stringr::str_remove(subsample_here, "^sample_size_"))
      )
  })
  
  if (nrow(pairwise_logfc_all) == 0) {
    message("  No usable log2FC values found for ", sample_type)
    next
  }
  
  log2fc_summary_df <- pairwise_logfc_all %>%
    group_by(sequencing_depth, tech1, tech2) %>%
    summarise(
      median_log2FC = median(log2FC, na.rm = TRUE),
      q1_log2FC = quantile(log2FC, 0.25, na.rm = TRUE),
      q3_log2FC = quantile(log2FC, 0.75, na.rm = TRUE),
      n_genes = dplyr::n(),
      .groups = "drop"
    )
  
  all_depths <- sort(unique(log2fc_summary_df$sequencing_depth))
  
  if (length(all_depths) == 0) {
    message("  No sequencing depths found for ", sample_type)
    next
  }
  
  all_pairs <- expand.grid(
    tech1 = technologies,
    tech2 = technologies,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    filter(match(tech1, technologies) < match(tech2, technologies))
  
  complete_grid <- tidyr::crossing(
    sequencing_depth = all_depths,
    all_pairs
  )
  
  plot_df_lower <- complete_grid %>%
    left_join(log2fc_summary_df, by = c("sequencing_depth", "tech1", "tech2")) %>%
    mutate(
      tech1 = factor(tech1, levels = technologies),
      tech2 = factor(tech2, levels = technologies)
    )
  
  base_depth_theme <- theme_bw() +
    theme(
      panel.grid.major.x = element_line(colour = "grey80", linewidth = 0.15),
      panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.15),
      panel.grid.minor = element_blank(),
      text = element_text(color = "black", size = 6),
      axis.text = element_text(size = 7, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.08, "cm")
    )
  
  make_blank_matrix_panel <- function(all_depths) {
    xlims <- range(all_depths, na.rm = TRUE)
    if (!all(is.finite(xlims))) xlims <- c(0, 1)
    
    ggplot() +
      scale_x_continuous(
        limits = xlims,
        labels = scales::label_number(scale = 1e-6, big.mark = "", suffix = "M"),
        expand = expansion(mult = c(0.02, 0.04))
      ) +
      scale_y_continuous(
        limits = c(-7, 7),
        breaks = seq(-8, 8, by = 2),
        expand = c(0, 0)
      ) +
      base_depth_theme +
      labs(x = NULL, y = NULL) +
      theme(
        panel.background = element_rect(fill = "grey90", color = NA)
      )
  }
  
  make_depth_panel <- function(panel_df, all_depths) {
    if (nrow(panel_df) == 0 || all(is.na(panel_df$median_log2FC))) {
      return(make_blank_matrix_panel(all_depths))
    }
    
    ggplot(panel_df, aes(x = sequencing_depth, y = median_log2FC, group = 1)) +
      geom_hline(
        yintercept = 0,
        linetype = "longdash",
        color = "black",
        linewidth = 0.5
      ) +
      geom_ribbon(
        aes(ymin = q1_log2FC, ymax = q3_log2FC),
        fill = "grey60",
        alpha = 0.5,
        na.rm = TRUE
      ) +
      geom_line(
        aes(y = q1_log2FC),
        linetype = "dotted",
        color = "black",
        linewidth = 0.3,
        na.rm = TRUE
      ) +
      geom_line(
        aes(y = q3_log2FC),
        linetype = "dotted",
        color = "black",
        linewidth = 0.3,
        na.rm = TRUE
      ) +
      geom_line(
        color = "black",
        linewidth = 0.5,
        na.rm = TRUE
      ) +
      geom_point(
        color = "black",
        shape = 18,
        size = 0.35,
        na.rm = TRUE
      ) +
      scale_x_continuous(
        labels = scales::label_number(scale = 1e-6, big.mark = "", suffix = "M"),
        expand = expansion(mult = c(0.02, 0.04))
      ) +
      scale_y_continuous(
        limits = c(-7, 7),
        breaks = seq(-8, 8, by = 2),
        expand = c(0, 0)
      ) +
      base_depth_theme +
      labs(x = NULL, y = NULL)
  }
  
  ntech <- length(technologies)
  plot_list <- vector("list", ntech * ntech)
  
  for (row in seq_len(ntech)) {
    for (col in seq_len(ntech)) {
      idx <- (row - 1) * ntech + col
      
      row_tech <- technologies[row]
      col_tech <- technologies[col]
      
      if (row > col) {
        panel_df <- plot_df_lower %>%
          filter(
            as.character(tech1) == col_tech,
            as.character(tech2) == row_tech
          )
        
        plot_list[[idx]] <- make_depth_panel(panel_df, all_depths)
        
      } else if (row < col) {
        panel_df <- plot_df_lower %>%
          filter(
            as.character(tech1) == row_tech,
            as.character(tech2) == col_tech
          ) %>%
          mutate(
            median_log2FC_orig = median_log2FC,
            q1_log2FC_orig = q1_log2FC,
            q3_log2FC_orig = q3_log2FC,
            median_log2FC = -median_log2FC_orig,
            q1_log2FC = -q3_log2FC_orig,
            q3_log2FC = -q1_log2FC_orig
          ) %>%
          select(-median_log2FC_orig, -q1_log2FC_orig, -q3_log2FC_orig)
        
        plot_list[[idx]] <- make_depth_panel(panel_df, all_depths)
        
      } else {
        plot_list[[idx]] <- make_blank_matrix_panel(all_depths)
      }
    }
  }
  
  matrix_plot <- GGally::ggmatrix(
    plots = plot_list,
    nrow = ntech,
    ncol = ntech,
    xAxisLabels = vapply(technologies, wrap_tech_label, character(1), width = 10),
    yAxisLabels = vapply(technologies, wrap_tech_label, character(1), width = 10)
  ) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
      strip.text = element_text(size = 7, lineheight = 0.9)
    )
  
  ggsave(
    file.path(outdir, paste0(sample_type, "_log2FC_sequencing_depth.png")),
    matrix_plot,
    width = 8,
    height = 8,
    dpi = 300
  )
}