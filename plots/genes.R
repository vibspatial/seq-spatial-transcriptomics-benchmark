# Libraries ---------------------------------------------------------------
library(tidyverse)
library(UpSetR)
library(scales)

# Root directories --------------------------------------------------------
root_dir <- "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed"
out_root <- file.path(root_dir, "plots/genes")

# Define look-up table ----------------------------------------------------
sample_metadata <- tribble(
  ~Sample,      ~sample_type,    ~Technology,       ~replicate,    ~file_path,
  "SCA013",     "Mouse lung",     "Visium",         NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/harpy/gene_total_counts.csv",
  "SCA014",     "Mouse liver",    "Visium",         NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/harpy/gene_total_counts.csv",
  "SCA016",     "Mouse brain",    "Visium",         NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA016/harpy/gene_total_counts.csv",
  "SCA015",     "Human lung",     "Visium",         NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/harpy/gene_total_counts.csv",
  
  "SCA001",     "Mouse lung",     "Visium_HD_FFPE", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA001/harpy/gene_total_counts.csv",
  "SCA002",     "Mouse liver",    "Visium_HD_FFPE", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA002/harpy/gene_total_counts.csv",
  "SCA003",     "Mouse brain",    "Visium_HD_FFPE", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA003/harpy/gene_total_counts.csv",
  "SCA010",     "Human lung",     "Visium_HD_FFPE", NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010/harpy/gene_total_counts.csv",
  
  "SCA017",     "Mouse lung",     "Visium_HD_FF",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA017/harpy/gene_total_counts.csv",
  "SCA018",     "Mouse liver",    "Visium_HD_FF",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA018/harpy/gene_total_counts.csv",
  "SCA019",     "Mouse brain",    "Visium_HD_FF",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA019/harpy/gene_total_counts.csv",
  "SCA020",     "Human lung",     "Visium_HD_FF",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA020/harpy/gene_total_counts.csv",
  
  "SCA021",     "Mouse liver",    "Visium_HD_3",    1,              "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA021/harpy/gene_total_counts.csv",
  "SCA023",     "Mouse liver",    "Visium_HD_3",    2,              "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA023/harpy/gene_total_counts.csv",
  "SCA024",     "Mouse brain",    "Visium_HD_3",    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA024/harpy/gene_total_counts.csv",
  "FPE010",     "Rice root",      "Visium_HD_3",    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/FPE010/harpy/gene_total_counts.csv",
  
  "Exp74",      "Mouse brain",    "BMKMANU_S1000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/exp74/harpy/gene_total_counts.csv",
  
  "Exp101",     "Mouse lung",     "BMKMANU_S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/harpy/gene_total_counts.csv",
  "Exp100",     "Mouse liver",    "BMKMANU_S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/harpy/gene_total_counts.csv",
  "Exp99",      "Mouse brain",    "BMKMANU_S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/harpy/gene_total_counts.csv",
  "Exp102",     "Human lung",     "BMKMANU_S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/harpy/gene_total_counts.csv",
  "psb",        "Rice root",      "BMKMANU_S3000",  NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/psb/harpy/gene_total_counts.csv",
  
  "SC12",       "Mouse lung",     "Curio_Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC12/harpy/gene_total_counts.csv",
  "SC13",       "Mouse liver",    "Curio_Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC13/harpy/gene_total_counts.csv",
  "SC11",       "Mouse brain",    "Curio_Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/harpy/gene_total_counts.csv",
  "SC14",       "Human lung",     "Curio_Seeker",   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC14/harpy/gene_total_counts.csv",
  
  "Exp93-sampleSPC002", "Mouse lung",  "Nova-ST",    NA, "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC002/harpy/gene_total_counts.csv",
  "Exp100-sampleSPC004","Mouse liver", "Nova-ST",    NA, "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp100-sampleSPC004/harpy/gene_total_counts.csv",
  "Exp65-sampleSPC014", "Mouse brain", "Nova-ST",    NA, "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp65-sampleSPC014/harpy/gene_total_counts.csv",
  "Exp93-sampleSPC022", "Human lung",  "Nova-ST",    NA, "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC022/harpy/gene_total_counts.csv",
  "Exp100-samplePSB",   "Rice root",   "Nova-ST",    NA, "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp100-samplePSB/harpy/gene_total_counts.csv",
  
  "spc002",     "Mouse lung",     "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc002/harpy/gene_total_counts.csv",
  "spc004",     "Mouse liver",    "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc004/harpy/gene_total_counts.csv",
  "spc014",     "Mouse brain",    "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc014/harpy/gene_total_counts.csv",
  "spc022",     "Human lung",     "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc022/harpy/gene_total_counts.csv",
  "psb_rice",   "Rice root",      "Stereo-Seq",     NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/psb/harpy/gene_total_counts.csv"
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

sample_type_levels <- c(
  "Rice root",
  "Human lung",
  "Mouse lung",
  "Mouse liver",
  "Mouse brain"
)


# Read all data -----------------------------------------------------------
read_gene_counts <- function(sample, path, sample_type, technology, replicate) {
  read.csv(path) %>%
    mutate(
      Sample = sample,
      sample_type = sample_type,
      Technology = technology,
      Replicate = replicate
    )
}

all_data <- pmap_dfr(
  sample_metadata,
  function(Sample, sample_type, Technology, replicate, file_path) {
    read_gene_counts(Sample, file_path, sample_type, Technology, replicate)
  }
) %>%
  mutate(
    Technology = recode(
      Technology,
      "Visium_HD_3" = "Visium HD 3'",
      "Visium_HD_FF" = "Visium HD (FF)",
      "Visium_HD_FFPE" = "Visium HD (FFPE)",
      "Visium" = "Visium (FFPE)",
      "Stereo-Seq" = "Stereo-Seq",
      "Nova-ST" = "Nova-ST",
      "Curio_Seeker" = "Curio Seeker",
      "BMKMANU_S3000" = "BMKMANU S3000",
      "BMKMANU_S1000" = "BMKMANU S1000"
    ),
    Technology = factor(Technology, levels = tech_levels),
    sample_type = factor(sample_type, levels = sample_type_levels)
  )

make_tech_label <- function(technology, replicate) {
  ifelse(
    is.na(replicate),
    as.character(technology),
    paste0(as.character(technology), "_", replicate)
  )
}

# UpSet prep --------------------------------------------------------------
prepare_presence_long <- function(sample_metadata_subset, all_data) {
  all_data %>%
    semi_join(sample_metadata_subset %>% select(Sample), by = "Sample") %>%
    transmute(
      gene,
      Technology,
      Replicate,
      tech_label = make_tech_label(Technology, Replicate),
      present = as.integer(total_counts > 0)
    )
}

make_presence_matrix <- function(df_long) {
  df_long %>%
    distinct(gene, tech_label, present) %>%
    pivot_wider(
      names_from = tech_label,
      values_from = present,
      values_fill = 0
    ) %>%
    arrange(gene)
}

plot_upset <- function(df_wide, sample_type, outdir) {
  library(dplyr)
  library(tibble)
  library(UpSetR)
  library(stringr)
  
  set_cols <- setdiff(colnames(df_wide), "gene")
  
  set_order <- unlist(lapply(tech_levels, function(tech) {
    set_cols[set_cols == tech | startsWith(set_cols, paste0(tech, "_"))]
  }))
  
  set_order <- unique(set_order)
  
  upset_data <- df_wide %>%
    column_to_rownames("gene") %>%
    as.data.frame()
  
  upset_data <- upset_data[, set_order, drop = FALSE]
  
  pretty_names <- colnames(upset_data) %>%
    str_replace("_(\\d+)$", " (\\1)")
  
  colnames(upset_data) <- pretty_names
  
  outfile <- file.path(outdir, paste0(sample_type, "_upset.png"))
  
  png(outfile, width = 2800, height = 1800, res = 300)
  
  par(
    mar = c(10, 12, 6, 4),
    mgp = c(3, 1, 0)
  )
  
  p <- upset(
    upset_data,
    sets = pretty_names,
    keep.order = TRUE,
    nsets = length(pretty_names),
    nintersects = 15,
    order.by = "freq",
    
    mainbar.y.label = "Gene intersections",
    sets.x.label = "Genes per sample",
    
    text.scale = c(
      2.0,  # intersection size title
      1.8,  # intersection size ticks
      1.7,  # set size title
      1.6,  # set size ticks
      1.5,  # set names
      1.5   # numbers above bars
    ),
    
    mb.ratio = c(0.65, 0.35),
    
    number.angles = 0,
    
    main.bar.color = "black",
    sets.bar.color = "black",
    matrix.color = "black",
    shade.color = "grey90",
    
    point.size = 4,
    line.size = 1.2
  )
  
  print(p)
  dev.off()
}

# Cumulative contribution -------------------------------------------------
get_common_genes <- function(df_long) {
  required_sets <- df_long %>%
    distinct(tech_label) %>%
    nrow()
  
  df_long %>%
    filter(present == 1) %>%
    distinct(gene, tech_label) %>%
    count(gene, name = "n_sets") %>%
    filter(n_sets == required_sets) %>%
    pull(gene)
}

compute_cumulative <- function(df_subset, gene_subset = NULL) {
  df_plot <- df_subset %>%
    filter(!is.na(normalized_counts), normalized_counts > 0) %>%
    mutate(
      tech_label = make_tech_label(Technology, Replicate)
    )
  
  if (!is.null(gene_subset)) {
    df_plot <- df_plot %>%
      filter(gene %in% gene_subset)
  }
  
  df_plot %>%
    group_by(tech_label, Technology, Replicate) %>%
    arrange(desc(normalized_counts), .by_group = TRUE) %>%
    mutate(
      rank = row_number(),
      cumulative = cumsum(normalized_counts)
    ) %>%
    ungroup()
}

plot_cumulative <- function(df_cum, sample_type, outdir, tech_colors,
                            log_scale = FALSE, y_label = "Cumalative contribution (%)", prefix = "cumulative") {
  df_plot <- df_cum %>%
    mutate(Technology = factor(Technology, levels = names(tech_colors)))
  
  label_df <- df_plot %>%
    filter(!is.na(Replicate)) %>%
    group_by(tech_label) %>%
    filter(rank == min(rank, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(Technology, Replicate) %>%
    mutate(
      y_label = (cumulative * 100) + (Replicate * 1.5)
    )
  
  p <- ggplot(
    df_plot,
    aes(
      x = rank,
      y = cumulative * 100,
      color = Technology,
      group = tech_label
    )
  ) +
    geom_line(size = 1) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      text = element_text(color = "black", size = 12),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(size = 0.2, color = "#DFDFDF"),
      panel.grid.major.y = element_line(size = 0.2, color = "#DFDFDF"),
      axis.text = element_text(color = "black", size = 11),
      axis.title.x = element_text(size = 12, margin = margin(r = 5)),
      axis.title.y = element_text(size = 11),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    labs(
      title = NULL,
      x = str_wrap("Number of genes", width = 25),
      y = str_wrap(y_label, width = 40),
      color = "Technology") +
    scale_color_manual(values = tech_colors) +
    guides(color = guide_legend(reverse = FALSE)) +
    geom_text(
      data = label_df,
      aes(y = y_label, label = Replicate),
      hjust = 1.2,
      vjust = 1.3,
      size = 2.5,
      show.legend = FALSE,
      color = "black"
    )
  
  if (log_scale) {
    p <- p +
      scale_x_log10(
        labels = scales::comma_format()
      )
  } else {
    p <- p +
      scale_x_continuous(
        labels = scales::comma_format(),
        expand = expansion(mult = c(0, 0.05))
      )
  }

  p <- p +
    scale_y_continuous(
      breaks = seq(0, 100, by = 10),
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.02))
    )
  
  suffix <- ifelse(log_scale, "log", "linear")
  outfile <- file.path(outdir, paste0(sample_type, "_", prefix, "_", suffix, ".png"))
  ggsave(outfile, p, width = 6, height = 4, dpi = 300)
}

get_gene_sets_list <- function(df_long) {
  df_long %>%
    filter(present == 1) %>%
    distinct(tech_label, gene) %>%
    group_by(tech_label) %>%
    summarise(genes = list(sort(gene)), .groups = "drop")
}

compute_pairwise_set_similarity <- function(df_long) {
  gene_sets <- get_gene_sets_list(df_long)
  set_names <- gene_sets$tech_label
  
  combn(set_names, 2, simplify = FALSE) %>%
    purrr::map_dfr(function(pair) {
      a <- gene_sets$genes[[match(pair[1], gene_sets$tech_label)]]
      b <- gene_sets$genes[[match(pair[2], gene_sets$tech_label)]]
      
      n_intersection <- length(intersect(a, b))
      n_union <- length(union(a, b))
      n_a <- length(a)
      n_b <- length(b)
      
      tibble(
        set_1 = pair[1],
        set_2 = pair[2],
        n_set_1 = n_a,
        n_set_2 = n_b,
        intersection = n_intersection,
        union = n_union,
        jaccard = n_intersection / n_union,
        overlap_coef = n_intersection / min(n_a, n_b)
      )
    })
}

plot_similarity_heatmap <- function(pairwise_similarity, sample_type, outdir,
                                    metric = c("jaccard", "overlap_coef"),
                                    suffix = NULL,
                                    legend_title = NULL) {
  library(viridis)
  library(tidyverse)
  
  metric <- match.arg(metric)
  
  if (is.null(suffix)) {
    suffix <- metric
  }
  
  if (is.null(legend_title)) {
    legend_title <- dplyr::case_when(
      metric == "jaccard" ~ "Jaccard",
      metric == "overlap_coef" ~ "Overlap",
      TRUE ~ metric
    )
  }
  
  all_sets <- sort(unique(c(pairwise_similarity$set_1, pairwise_similarity$set_2)))
  
  metric_sym <- rlang::sym(metric)
  metric_rev <- paste0(metric, "_rev")
  
  heatmap_df <- expand.grid(set_1 = all_sets, set_2 = all_sets, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    left_join(
      pairwise_similarity %>% select(set_1, set_2, !!metric_sym),
      by = c("set_1", "set_2")
    ) %>%
    left_join(
      pairwise_similarity %>%
        select(set_1, set_2, !!metric_sym) %>%
        rename(
          set_1 = set_2,
          set_2 = set_1,
          !!metric_rev := !!metric_sym
        ),
      by = c("set_1", "set_2")
    ) %>%
    mutate(
      value = coalesce(.data[[metric]], .data[[metric_rev]]),
      value = ifelse(set_1 == set_2, 1, value)
    ) %>%
    select(set_1, set_2, value) %>%
    mutate(
      set_1 = factor(set_1, levels = all_sets),
      set_2 = factor(set_2, levels = rev(all_sets))
    )
  
  p <- ggplot(heatmap_df, aes(x = set_1, y = set_2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
    scale_fill_viridis(
      option = "viridis",
      limits = c(0, 1),
      name = legend_title
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(color = "black", size = 12),
      axis.text = element_text(color = "black", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_text(size = 12, margin = margin(r = 5)),
      axis.title.y = element_text(size = 12),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = NULL
    )
  
  ggsave(
    file.path(outdir, paste0(sample_type, "_", suffix, "_heatmap.png")),
    p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

make_common_gene_expression_table <- function(df_subset, df_long) {
  common_genes <- get_common_genes(df_long)
  
  df_common <- df_subset %>%
    filter(gene %in% common_genes, !is.na(normalized_counts)) %>%
    mutate(
      tech_label = make_tech_label(Technology, Replicate)
    ) %>%
    group_by(gene, tech_label) %>%
    summarise(
      normalized_counts = sum(normalized_counts, na.rm = TRUE),
      .groups = "drop"
    )  %>%
    group_by(tech_label) %>%
    mutate(normalized_counts = normalized_counts / sum(normalized_counts, na.rm = TRUE)) %>%
    ungroup()
  
  df_wide <- df_common %>%
    pivot_wider(
      names_from = tech_label,
      values_from = normalized_counts,
      values_fill = 0
    )
  
  sample_cols <- setdiff(colnames(df_wide), "gene")
  
  df_wide <- df_wide %>%
    mutate(
      average_normalized_counts = rowMeans(across(all_of(sample_cols)), na.rm = TRUE)
    )
  
  df_wide <- df_wide %>%
    mutate(
      expression_pentile = dplyr::ntile(average_normalized_counts, 5),
      expression_pentile = recode(
        expression_pentile,
        `1` = "very low",
        `2` = "low",
        `3` = "mid",
        `4` = "high",
        `5` = "very high"
      )
    ) %>%
    arrange(desc(average_normalized_counts))
  
  df_wide
}

make_common_gene_sample_sums <- function(df_subset, df_long) {
  common_genes <- get_common_genes(df_long)
  
  df_subset %>%
    filter(gene %in% common_genes, !is.na(normalized_counts)) %>%
    mutate(
      tech_label = make_tech_label(Technology, Replicate)
    ) %>%
    group_by(tech_label, Technology, Replicate) %>%
    summarise(
      common_genes_sum_normalized_counts = sum(normalized_counts, na.rm = TRUE),
      n_common_genes_detected = sum(normalized_counts > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(factor(Technology, levels = tech_levels), Replicate)
}

plot_common_gene_sums_bar <- function(common_gene_sums, sample_type, outdir, tech_colors) {
  plot_order <- common_gene_sums %>%
    mutate(
      Technology = factor(Technology, levels = names(tech_colors)),
      Replicate_order = ifelse(is.na(Replicate), 0, Replicate)
    ) %>%
    arrange(Technology, Replicate_order) %>%
    pull(tech_label) %>%
    unique()
  
  df_plot <- common_gene_sums %>%
    mutate(
      Technology = factor(Technology, levels = names(tech_colors)),
      Replicate_order = ifelse(is.na(Replicate), 0, Replicate),
      Sample = factor(tech_label, levels = rev(plot_order)),
      Sample_num = as.numeric(Sample),
      y_val = common_genes_sum_normalized_counts * 100
    )
  
  dodge_width <- 0.9
  max_val <- max(df_plot$y_val, na.rm = TRUE)
  
  p <- ggplot(df_plot, aes(x = Sample, y = y_val, fill = Technology)) +
    geom_bar(
      aes(group = Technology),
      stat = "identity",
      position = position_dodge2(width = dodge_width, preserve = "single", padding = 0),
      color = "black",
      width = 0.9
    ) +
    geom_text(
      aes(label = sprintf("%.1f", y_val), group = Technology),
      position = position_dodge2(width = dodge_width, preserve = "single"),
      hjust = -0.1,
      vjust = 0.5,
      size = 2.6
    ) +
    geom_text(
      data = df_plot %>% filter(!is.na(Replicate)),
      aes(label = Replicate, group = Technology),
      position = position_dodge2(width = dodge_width, preserve = "single"),
      y = 0.005 * max_val,
      hjust = 0,
      vjust = 0.5,
      size = 2.6
    ) +
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
    scale_x_discrete(labels = function(x) {
      as.character(df_plot$Technology[match(x, df_plot$Sample)])
    }) +
    labs(
      title = NULL,
      y = str_wrap("Contribution to unique UMI counts in tissue by common genes (%)", width = 35),
      x = NULL,
      fill = "Technology"
    ) +
    scale_fill_manual(values = tech_colors) +
    guides(fill = guide_legend(reverse = FALSE)) +
    scale_y_continuous(
      breaks = seq(0, 101, by = 10),
      expand = expansion(mult = c(0, 0.1))
    )
  
  n_bars <- length(unique(df_plot$Sample))
  plot_height <- ( 1 + n_bars) * 0.5
  
  ggsave(
    file.path(outdir, paste0(sample_type, "_common_gene_sum_normalized_counts_barplot.png")),
    plot = p,
    width = 8,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )
}
plot_gene_counts_bar <- function(df_cum, sample_type, outdir, tech_colors, prefix = "gene_counts") {
  
  df_plot <- df_cum %>%
    group_by(tech_label, Technology, Replicate) %>%
    summarise(
      n_genes = max(rank, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Technology = factor(Technology, levels = names(tech_colors)),
      Replicate_order = ifelse(is.na(Replicate), 0, Replicate)
    ) %>%
    arrange(Technology, Replicate_order) %>%
    mutate(
      Sample = factor(tech_label, levels = rev(unique(tech_label)))
    )
  
  max_val <- max(df_plot$n_genes, na.rm = TRUE)
  
  p <- ggplot(df_plot, aes(x = Sample, y = n_genes, fill = Technology)) +

    geom_col(color = "black", width = 0.9) +
    
    # gene count labels
    geom_text(
      aes(label = scales::comma(n_genes)),
      hjust = -0.1,
      vjust = 0.5,
      size = 2.8
    ) +
    
    # replicate labels
    geom_text(
      data = df_plot %>% filter(!is.na(Replicate)),
      aes(label = Replicate),
      y = 0.03 * max_val,
      hjust = 0,
      vjust = 0.5,
      size = 3.2
    ) +
    
    coord_flip(clip = "off") +
    
    scale_fill_manual(values = tech_colors) +
    
    scale_y_continuous(
      limits = c(0, 40000),
    ) +
    
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 15, 0, 0)
    )
  
  outfile <- file.path(outdir, paste0(sample_type, "_", prefix, ".png"))
  
  ggsave(
    outfile,
    p,
    width = 2.2,
    height = 0.18 * nrow(df_plot),
    dpi = 300
  )
}

# Main loop ---------------------------------------------------------------
for (current_sample_type in levels(all_data$sample_type)) {
  message("Processing ", current_sample_type)
  
  sample_metadata_subset <- sample_metadata %>%
    filter(sample_type == current_sample_type)
  
  if (nrow(sample_metadata_subset) < 2) next
  
  outdir <- file.path(out_root, current_sample_type)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # UpSet
  df_long <- prepare_presence_long(sample_metadata_subset, all_data)
  presence_matrix <- make_presence_matrix(df_long)
  
  write.csv(
    presence_matrix,
    file.path(outdir, paste0(current_sample_type, "_presence_matrix.csv")),
    row.names = FALSE
  )
  
  plot_upset(presence_matrix, current_sample_type, outdir)
  
  # All genes cumulative contribution
  df_subset <- all_data %>%
    semi_join(sample_metadata_subset %>% select(Sample), by = "Sample")
  
  df_cum_all <- compute_cumulative(df_subset)
  
  plot_cumulative(
    df_cum_all, current_sample_type, outdir, tech_colors,
    log_scale = FALSE,
    prefix = "cumulative_all_genes",
    y_label = "Cumulative contribution to unique UMI counts in tissue (%)"
  )
  plot_cumulative(
    df_cum_all, current_sample_type, outdir, tech_colors,
    log_scale = TRUE,
    prefix = "cumulative_all_genes",
    y_label = "Cumulative contribution to unique UMI counts in tissue (%)"
  )
  
  plot_gene_counts_bar(
    df_cum_all,
    current_sample_type,
    outdir,
    tech_colors,
    prefix = "gene_counts_all"
  )
  
  # Common and non-common genes cumulative contribution
  common_genes <- get_common_genes(df_long)
  non_common_genes <- setdiff(unique(df_subset$gene), common_genes)
  
  write.csv(
    tibble(gene = sort(common_genes)),
    file.path(outdir, paste0(current_sample_type, "_common_genes.csv")),
    row.names = FALSE
  )
  
  write.csv(
    tibble(gene = sort(non_common_genes)),
    file.path(outdir, paste0(current_sample_type, "_non_common_genes.csv")),
    row.names = FALSE
  )
  
  df_cum_common <- compute_cumulative(df_subset, gene_subset = common_genes)
  df_cum_non_common <- compute_cumulative(df_subset, gene_subset = non_common_genes)
  
  plot_cumulative(
    df_cum_common, current_sample_type, outdir, tech_colors,
    log_scale = FALSE,
    prefix = "cumulative_common_genes",
    y_label = "Cumulative contribution to unique UMI counts in tissue by common genes (%)"
  )
  plot_cumulative(
    df_cum_common, current_sample_type, outdir, tech_colors,
    log_scale = TRUE,
    prefix = "cumulative_common_genes",
    y_label = "Cumulative contribution to unique UMI counts in tissue by common genes (%)"
  )
  
  plot_cumulative(
    df_cum_non_common, current_sample_type, outdir, tech_colors,
    log_scale = FALSE,
    prefix = "cumulative_non_common_genes",
    y_label = "Cumulative contribution to unique UMI counts in tissue by non-common genes (%)"
  )
  plot_cumulative(
    df_cum_non_common, current_sample_type, outdir, tech_colors,
    log_scale = TRUE,
    prefix = "cumulative_non_common_genes",
    y_label = "Cumulative contribution to unique UMI counts in tissue by non-common genes (%)"
  )
  
  plot_gene_counts_bar(df_cum_common, current_sample_type, outdir, tech_colors, "gene_counts_common")
  plot_gene_counts_bar(df_cum_non_common, current_sample_type, outdir, tech_colors, "gene_counts_non_common")
  
  # Similarity
  pairwise_similarity <- compute_pairwise_set_similarity(df_long)
  
  write.csv(
    pairwise_similarity,
    file.path(outdir, paste0(current_sample_type, "_pairwise_set_similarity.csv")),
    row.names = FALSE
  )
  
  plot_similarity_heatmap(
    pairwise_similarity,
    current_sample_type,
    outdir,
    metric = "jaccard",
    suffix = "jaccard",
    legend_title = "Jaccard"
  )
  
  plot_similarity_heatmap(
    pairwise_similarity,
    current_sample_type,
    outdir,
    metric = "overlap_coef",
    suffix = "overlap",
    legend_title = "Overlap"
  )
  
  # Common-gene expression table
  common_gene_expression <- make_common_gene_expression_table(df_subset, df_long)
  
  write.csv(
    common_gene_expression,
    file.path(outdir, paste0(current_sample_type, "_common_gene_expression_table.csv")),
    row.names = FALSE
  )
  
  # Per-sample sums over common genes
  common_gene_sums <- make_common_gene_sample_sums(df_subset, df_long)
  
  plot_common_gene_sums_bar(
    common_gene_sums,
    current_sample_type,
    outdir,
    tech_colors
  )
}
