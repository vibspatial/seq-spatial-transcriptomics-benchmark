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
  ~Sample,      ~sample_type,    ~Technology,       ~area,   ~replicate,     ~bam_dir,
  "SCA013",     "Mouse_lung",     "Visium",         15.46,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/bam_output",
  "SCA014",     "Mouse_liver",    "Visium",         27.97,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/bam_output",
  "SCA016",     "Mouse_brain",    "Visium",         24.18,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA016/bam_output",
  "SCA015",     "Human_lung",     "Visium",         34.37,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/bam_output",
  
  "SCA001",     "Mouse_lung",     "Visium_HD_FFPE", 15.36,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA001/bam_output",
  "SCA002",     "Mouse_liver",    "Visium_HD_FFPE", 29.54,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA002/bam_output",
  "SCA003",     "Mouse_brain",    "Visium_HD_FFPE", 24.84,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA003/bam_output",
  #"SCA005",     "Mouse_colon",    "Visium_HD_FFPE", 6.97,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA005/bam_output",
  "SCA010",     "Human_lung",     "Visium_HD_FFPE", 33.3,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010/bam_output",
  
  "SCA017",     "Mouse_lung",     "Visium_HD_FF",   10.4,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA017/bam_output",
  "SCA018",     "Mouse_liver",    "Visium_HD_FF",   21.19,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA018/bam_output",
  "SCA019",     "Mouse_brain",    "Visium_HD_FF",   19.04,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA019/bam_output",
  "SCA020",     "Human_lung",     "Visium_HD_FF",   40.05,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA020/bam_output",

  "SCA021",     "Mouse_liver",    "Visium_HD_3",    16.06,   1,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA021/bam_output",
  "SCA023",     "Mouse_liver",    "Visium_HD_3",    9.49,    2,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA023/bam_output",
  "SCA024",     "Mouse_brain",    "Visium_HD_3",    13.03,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA024/bam_output",
  "FPE010",     "Rice_root",      "Visium_HD_3",    0.48,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/FPE010/bam_output",
  
  "Exp74",      "Mouse_brain",    "BMKMANU_S1000",  21.8,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/exp74/bam_output",
  
  "Exp101",     "Mouse_lung",     "BMKMANU_S3000",  13.01,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/bam_output",
  "Exp100",     "Mouse_liver",    "BMKMANU_S3000",  19.6,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/bam_output",
  "Exp99",      "Mouse_brain",    "BMKMANU_S3000",  20.94,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/bam_output",
  "Exp102",     "Human_lung",     "BMKMANU_S3000",  25.65,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/bam_output",
  "psb",        "Rice_root",      "BMKMANU_S3000",  2.15,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/psb/bam_output",
  
  "SC12",       "Mouse_lung",     "Curio_Seeker",   8.57,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC12/bam_output",
  "SC13",       "Mouse_liver",    "Curio_Seeker",   8.57,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC13/bam_output",
  "SC11",       "Mouse_brain",    "Curio_Seeker",   8.61,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/bam_output",
  "SC14",       "Human_lung",     "Curio_Seeker",   8.53,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC14/bam_output",

  "Exp93-sampleSPC002",     "Mouse_lung",     "Nova-ST",    9.27,        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC002/bam_output",
  "Exp100-sampleSPC004",    "Mouse_liver",    "Nova-ST",    19.23,       NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp100-sampleSPC004/bam_output",
  "Exp65-sampleSPC014",     "Mouse_brain",    "Nova-ST",    20.27,       NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp65-sampleSPC014/bam_output",
  "Exp93-sampleSPC022",     "Human_lung",     "Nova-ST",    47.07,       NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC022/bam_output",
  "Exp100-samplePSB",       "Rice_root",      "Nova-ST",    0.45,        NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp100-samplePSB/bam_output",
  
  "spc002",     "Mouse_lung",     "Stereo-Seq",     13.47,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc002/bam_output",
  "spc004",     "Mouse_liver",    "Stereo-Seq",     17.21,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc004/bam_output",
  "spc014",     "Mouse_brain",    "Stereo-Seq",     22.09,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc014/bam_output",
  "spc022",     "Human_lung",     "Stereo-Seq",     44.24,   NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc022/bam_output",
  "psb_rice",   "Rice_root",      "Stereo-Seq",     0.98,    NA,             "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/psb/bam_output",
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
read_sample_metrics <- function(bam_dir, sample_name, sample_type, technology, area, replicate) {
  
  preseq_cols <- list(
    c_curve = c(
      "included_reads",
      "distinct",
      "corrected_reads"
    ),
    lc_extrap = c(
      "included_reads",
      "distinct",
      "lc_extrap_lower_ci",
      "lc_extrap_upper_ci",
      "corrected_reads"
    ),
    bound_pop = c(
      "bound_pop_estimation",
      "bound_pop_lower_ci",
      "bound_pop_upper_ci"
    )
  )
  
  read_preseq <- function(path, method) {
    stopifnot(method %in% names(preseq_cols))
    
    if (!file.exists(path)) {
      return(NULL)
    }
    
    cols <- preseq_cols[[method]]
    
    df <- read_tsv(
      path,
      col_names = cols,
      skip = 1,
      show_col_types = FALSE
    ) %>%
      tidyr::drop_na(all_of(cols))
    
    # Defensive check
    if (ncol(df) != length(cols)) {
      stop(
        "Column mismatch in ", path,
        " (expected ", length(cols),
        ", got ", ncol(df), ")"
      )
    }
    
    df
  }

  
  files <- tribble(
    ~path, ~method,      ~scope,
    file.path(bam_dir, "preseq", "c_curve_all.txt"),      "c_curve",   "all",
    file.path(bam_dir, "preseq", "c_curve_tissue.txt"),   "c_curve",   "tissue",
    file.path(bam_dir, "preseq", "lc_extrap_all.txt"),    "lc_extrap", "all",
    file.path(bam_dir, "preseq", "lc_extrap_tissue.txt"), "lc_extrap", "tissue",
    file.path(bam_dir, "preseq", "bound_pop_all.txt"),    "bound_pop", "all",
    file.path(bam_dir, "preseq", "bound_pop_tissue.txt"), "bound_pop", "tissue"
  )
  
  df <- files %>%
    mutate(
      data = purrr::map2(path, method, read_preseq)
    ) %>%
    dplyr::filter(!purrr::map_lgl(data, is.null)) %>%
    tidyr::unnest(data)
  
  if (nrow(df) == 0) {
    return(tibble::tibble())
  }
  
  df %>%
    mutate(
      Sample = sample_name,
      sample_type = sample_type,
      Technology = technology,
      bam_dir = bam_dir,
      area = area,
      replicate = replicate
    )
}

# Load all samples --------------------------------------------------------
all_data <- map_dfr(1:nrow(sample_metadata), function(i) {
  row <- sample_metadata[i, ]
  read_sample_metrics(row$bam_dir, row$Sample, row$sample_type, row$Technology, row$area, row$replicate)
})
all_data <- all_data %>%
  mutate(
    distinct_per_mm2 = distinct / area,
    lc_extrap_lower_ci_per_mm2 = lc_extrap_lower_ci / area,
    lc_extrap_upper_ci_per_mm2 = lc_extrap_upper_ci / area,
    bound_pop_estimation_per_mm2 = bound_pop_estimation / area,
    bound_pop_lower_ci_per_mm2 = bound_pop_lower_ci / area,
    bound_pop_upper_ci_per_mm2 = bound_pop_upper_ci / area,
    sequencing_saturation_reads =
      if_else(corrected_reads > 0,
              (1 - (distinct / corrected_reads)) * 100,
              0),
    sequencing_saturation_UMIs =
      if_else(included_reads > 0,
              (1 - (distinct / included_reads)) * 100,
              0)
  ) %>%
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
make_line_plot <- function(df, xcol, ycol, outpath, color_by, palette, title, xlab, ylab) {
  df_plot <- df %>%
    rename(x_val = !!sym(xcol)) %>%
    rename(y_val = !!sym(ycol))
  
  use_area_corrected <- ycol == "distinct_per_mm2"
  percentage <- grepl("percentage|sequencing_saturation", ycol)
  
  # Per-sample deterministic jitter for bound_pop / labels
  sample_jitter_df <- df_plot %>%
    distinct(Sample, replicate) %>%
    arrange(Sample) %>%
    mutate(
      jitter_offset = seq(-0.15e9, 0.15e9, length.out = n())
    )
  
  max_val <- max(df_plot$y_val, na.rm = TRUE)
  y_labels_func <- if (max_val > 1e6) {
    scales::label_number(scale = 1e-6, big.mark = "", suffix = "M")
  } else if (percentage) {
    scales::label_number(accuracy = 1)
  } else {
    scales::comma
  }
  
  ribbon_data <- NULL
  if (!percentage) {
    ribbon_data <- df_plot %>%
      filter(method == "lc_extrap", x_val < 4e9) %>%
      mutate(
        ymin = if (use_area_corrected) lc_extrap_lower_ci_per_mm2 else lc_extrap_lower_ci,
        ymax = if (use_area_corrected) lc_extrap_upper_ci_per_mm2 else lc_extrap_upper_ci
      ) %>%
      group_by(Sample) %>%
      arrange(x_val) %>%
      mutate(
        ymin_smooth = predict(loess(ymin ~ x_val, span = 0.3)),
        ymax_smooth = predict(loess(ymax ~ x_val, span = 0.3))
      ) %>%
      ungroup()
  }
  
  bound_pop_df <- NULL
  if (!percentage) {
    bound_pop_df <- df_plot %>%
      filter(method == "bound_pop") %>%
      distinct(
        Sample, Technology, sample_type, `Sample type`, replicate,
        bound_pop_estimation, bound_pop_lower_ci, bound_pop_upper_ci,
        bound_pop_estimation_per_mm2, bound_pop_lower_ci_per_mm2, bound_pop_upper_ci_per_mm2,
        .keep_all = TRUE
      ) %>%
      left_join(sample_jitter_df, by = c("Sample", "replicate")) %>%
      mutate(
        x_val = 4.2e9,
        x_bound = x_val + jitter_offset,
        y_val = if (use_area_corrected) bound_pop_estimation_per_mm2 else bound_pop_estimation,
        ymin = if (use_area_corrected) bound_pop_lower_ci_per_mm2 else bound_pop_lower_ci,
        ymax = if (use_area_corrected) bound_pop_upper_ci_per_mm2 else bound_pop_upper_ci
      )
  }

  # Compute y-axis upper limit from all plotted y-values
  y_max_plot <- if (percentage) {
    100
  } else {
    max_lc <- df_plot %>%
      filter(method == "lc_extrap", x_val < 4e9) %>%
      summarise(val = max(y_val, na.rm = TRUE)) %>%
      pull(val)
    
    max_bp <- bound_pop_df %>%
      summarise(val = max(y_val, na.rm = TRUE)) %>%
      pull(val)
    
    1.2 * max(max_lc, max_bp, na.rm = TRUE)
  }
  
  p <- ggplot(
    df_plot,
    aes(
      x = x_val,
      y = y_val,
      color = !!sym(color_by),
      group = interaction(Sample, method)
    )
  ) +
    (if (!percentage) {
      geom_ribbon(
        data = ribbon_data,
        aes(
          x = x_val,
          ymin = ymin_smooth,
          ymax = ymax_smooth,
          fill = !!sym(color_by),
          group = Sample
        ),
        alpha = 0.25,
        color = NA,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    } else NULL) +
    geom_line(
      data = df_plot %>% filter(method != "bound_pop") %>% filter(x_val < 4e9),
      aes(linetype = method),
      linewidth = 1.1
    ) +
    scale_linetype_manual(
      values = c(c_curve = "solid", lc_extrap = "longdash"),
      name = "Method"
    ) +
    (if (!percentage) {
      geom_vline(
        xintercept = 4e9,
        linetype = "dashed",
        color = "black",
        linewidth = 0.5
      )
    } else NULL) +
    (if (!percentage) {
      geom_vline(
        xintercept = 4.4e9,
        linetype = "dashed",
        color = "black",
        linewidth = 0.5
      )
    } else NULL) +
    (if (!percentage) {
      geom_errorbar(
        data = bound_pop_df,
        aes(
          x = x_bound,
          ymin = ymin,
          ymax = ymax,
          color = !!sym(color_by)
        ),
        width = 0.04e9,
        linewidth = 0.6,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    } else NULL) +
    (if (!percentage) {
      geom_point(
        data = bound_pop_df,
        aes(
          x = x_bound,
          y = y_val,
          color = !!sym(color_by)
        ),
        shape = 18,
        size = 2.4,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    } else NULL) +
    geom_point(
      data = df_plot %>%
        filter(method == "c_curve") %>%
        group_by(Sample) %>%
        slice_max(order_by = x_val, n = 1, with_ties = FALSE) %>%
        ungroup(),
      aes(x = x_val, y = y_val, color = !!sym(color_by)),
      shape = 21,
      size = 1.6,
      fill = "white",
      stroke = 1.2
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      text = element_text(color = "black", size = 12),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.2, color = "#DFDFDF"),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "#DFDFDF"),
      axis.text = element_text(color = "black", size = 11),
      axis.title.x = element_text(size = 12, margin = margin(r = 5)),
      axis.title.y = element_text(size = 12),
      plot.margin = margin(10, 40, 10, 10)
    ) +
    labs(title = NULL, x = xlab, y = str_wrap(ylab, width = 35), color = color_by) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, big.mark = "", suffix = "M"),
      breaks = pretty_breaks(n = 5),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.08))
    ) +
    scale_y_continuous(
      breaks = if (percentage) seq(0, 100, by = 10) else pretty_breaks(n = 6),
      labels = y_labels_func,
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(
      ylim = c(0, y_max_plot),
      clip = "off"
    ) +
    guides(
      color = guide_legend(reverse = TRUE, order = 2),
      linetype = guide_legend(override.aes = list(linewidth = 0.55), order = 1)
    )
  
  ggsave(outpath, plot = p, width = 8, height = 4)
}

# Generate plots per sample -----------------------------------------------
samples <- distinct(all_data, Sample, bam_dir)
for (i in seq_len(nrow(samples))) {
  sub <- filter(all_data, Sample == samples$Sample[i])
  outdir <- file.path(samples$bam_dir[i], "plots", "preseq")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  sub_all <- filter(sub, scope == "all")
  make_line_plot(sub_all, "corrected_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_all_reads.png"), "Technology", tech_colors, "Unique UMI counts", "Reads", "Unique UMI counts")
  make_line_plot(sub_all, "included_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_UMI_reads.png"), "Technology", tech_colors, "Unique UMI counts", "UMI counts", "Unique UMI counts")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_reads", file.path(outdir, "preseq_sequencing_saturation_reads_vs_all_reads.png"), "Technology", tech_colors, "Sequencing saturation (%)", "Reads", "Sequencing saturation (%)")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_all_reads.png"), "Technology", tech_colors, "Sequencing saturation (%)", "Reads", "Sequencing saturation (%)")
  make_line_plot(sub_all, "included_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_UMI_reads.png"), "Technology", tech_colors, "Sequencing saturation (%)", "UMI counts", "Sequencing saturation (%)")
  
  sub_tissue <- filter(sub, scope == "tissue")
  make_line_plot(sub_tissue, "corrected_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_all_reads.png"), "Technology", tech_colors, "Unique UMI counts detected in tissue", "Reads", "Unique UMI counts detected in tissue")
  make_line_plot(sub_tissue, "included_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_UMI_reads.png"), "Technology", tech_colors, "Unique UMI counts detected in tissue", "UMI counts", "Unique UMI counts detected in tissue")
  make_line_plot(sub_tissue, "corrected_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_all_reads.png"), "Technology", tech_colors, "Unique UMI counts detected in tissue/mm²", "Reads", "Unique UMI counts detected in tissue/mm²")
  make_line_plot(sub_tissue, "included_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_UMI_reads.png"), "Technology", tech_colors, "Unique UMI counts detected in tissue/mm²", "UMI counts", "Unique UMI counts detected in tissue/mm²")
}

# Generate plots per technology -------------------------------------------
tech_list <- unique(all_data$Technology)
for (tech in tech_list) {
  sub <- all_data[all_data$Technology == tech, ]
  
  outdir <- file.path(tech_outdirs[[tech]], "preseq")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  sub_all <- filter(sub, scope == "all")
  make_line_plot(sub_all, "corrected_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_all_reads.png"), "Sample type", sample_type_colors, "Unique UMI counts", "Reads", "Unique UMI counts")
  make_line_plot(sub_all, "included_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_UMI_reads.png"), "Sample type", sample_type_colors, "Unique UMI counts", "UMI counts", "Unique UMI counts")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_reads", file.path(outdir, "preseq_sequencing_saturation_reads_vs_all_reads.png"), "Sample type", sample_type_colors, "Sequencing saturation (%)", "Reads", "Sequencing saturation (%)")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_all_reads.png"), "Sample type", sample_type_colors, "Sequencing saturation (%)", "Reads", "Sequencing saturation (%)")
  make_line_plot(sub_all, "included_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_UMI_reads.png"), "Sample type", sample_type_colors, "Sequencing saturation (%)", "UMI counts", "Sequencing saturation (%)")
  
  sub_tissue <- filter(sub, scope == "tissue")
  make_line_plot(sub_tissue, "corrected_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_all_reads.png"), "Sample type", sample_type_colors, "Unique UMI counts detected in tissue", "Reads", "Unique UMI counts detected in tissue")
  make_line_plot(sub_tissue, "included_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_UMI_reads.png"), "Sample type", sample_type_colors, "Unique UMI counts detected in tissue", "UMI counts", "Unique UMI counts detected in tissue")
  make_line_plot(sub_tissue, "corrected_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_all_reads.png"), "Sample type", sample_type_colors, "Unique UMI counts detected in tissue/mm²", "Reads", "Unique UMI counts detected in tissue/mm²")
  make_line_plot(sub_tissue, "included_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_UMI_reads.png"), "Sample type", sample_type_colors, "Unique UMI counts detected in tissue/mm²", "UMI counts", "Unique UMI counts detected in tissue/mm²")
}

# Generate plots per sample type ------------------------------------------
sample_type_list <- unique(all_data$`Sample type`)
for (sample_type in sample_type_list) {
  sub <- all_data[all_data$`Sample type` == sample_type, ]
  
  outdir <- file.path("/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/plots/bam_plots", sample_type, "preseq")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  sub_all <- filter(sub, scope == "all")
  make_line_plot(sub_all, "corrected_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_all_reads.png"), "Technology", tech_colors, "Unique UMI counts", "Reads", "Unique UMI counts")
  make_line_plot(sub_all, "included_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_UMI_reads.png"), "Technology", tech_colors, "Unique UMI counts", "UMI counts", "Unique UMI counts")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_reads", file.path(outdir, "preseq_sequencing_saturation_reads_vs_all_reads.png"), "Technology", tech_colors, "Sequencing saturation (%)", "Reads", "Sequencing saturation (%)")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_all_reads.png"), "Technology", tech_colors, "Sequencing saturation (%)", "Reads", "Sequencing saturation (%)")
  make_line_plot(sub_all, "included_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_UMI_reads.png"), "Technology", tech_colors, "Sequencing saturation (%)", "UMI counts", "Sequencing saturation (%)")
  
  sub_tissue <- filter(sub, scope == "tissue")
  make_line_plot(sub_tissue, "corrected_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_all_reads.png"), "Technology", tech_colors, "Unique UMI counts detected in tissue", "Reads", "Unique UMI counts detected in tissue")
  make_line_plot(sub_tissue, "included_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_UMI_reads.png"), "Technology", tech_colors, "Unique UMI counts detected in tissue", "UMI counts", "Unique UMI counts detected in tissue")
  make_line_plot(sub_tissue, "corrected_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_all_reads.png"), "Technology", tech_colors, "Unique UMI counts detected in tissue/mm²", "Reads", "Unique UMI counts detected in tissue/mm²")
  make_line_plot(sub_tissue, "included_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_UMI_reads.png"), "Technology", tech_colors, "Unique UMI counts detected in tissue/mm²", "UMI counts", "Unique UMI counts detected in tissue/mm²")
}
