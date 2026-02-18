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
  ~Sample,      ~sample_type,    ~Technology,       ~area,    ~bam_dir,
  "SCA013",     "Mouse_lung",     "Visium",         15.46,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA013/bam_output",
  "SCA014",     "Mouse_liver",    "Visium",         27.97,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA014/bam_output",
  "SCA016",     "Mouse_brain",    "Visium",         24.18,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA016/bam_output",
  "SCA015",     "Human_lung",     "Visium",         34.37,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/bam_output",
  
  "SCA001",     "Mouse_lung",     "Visium_HD_FFPE", 15.36,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA001/bam_output",
  "SCA002",     "Mouse_liver",    "Visium_HD_FFPE", 29.54,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA002/bam_output",
  "SCA003",     "Mouse_brain",    "Visium_HD_FFPE", 24.84,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA003/bam_output",
  "SCA005",     "Mouse_colon",    "Visium_HD_FFPE", 6.97,    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA005/bam_output",
  "SCA010",     "Human_lung",     "Visium_HD_FFPE", 33.3,    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010/bam_output",
  
  "SCA017",     "Mouse_lung",     "Visium_HD_FF",   10.4,    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA017/bam_output",
  "SCA018",     "Mouse_liver",    "Visium_HD_FF",   21.19,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA018/bam_output",
  "SCA019",     "Mouse_brain",    "Visium_HD_FF",   19.04,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA019/bam_output",
  "SCA020",     "Human_lung",     "Visium_HD_FF",   40.05,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA020/bam_output",
  
  "Exp74",      "Mouse_brain",    "BMKMANU_S1000",  21.8,    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/exp74/bam_output",
  
  "Exp101",     "Mouse_lung",     "BMKMANU_S3000",  13.01,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp101/bam_output",
  "Exp100",     "Mouse_liver",    "BMKMANU_S3000",  19.6,    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp100/bam_output",
  "Exp99",      "Mouse_brain",    "BMKMANU_S3000",  20.94,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp99/bam_output",
  "Exp102",     "Human_lung",     "BMKMANU_S3000",  25.65,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/Exp102/bam_output",
  
  "SC12",       "Mouse_lung",     "Curio_Seeker",   8.57,    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC12/bam_output",
  "SC13",       "Mouse_liver",    "Curio_Seeker",   8.57,    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC13/bam_output",
  "SC11",       "Mouse_brain",    "Curio_Seeker",   8.61,    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC11/bam_output",
  "SC14",       "Human_lung",     "Curio_Seeker",   8.53,    "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC14/bam_output",

  "Exp93-sampleSPC002",     "Mouse_lung",     "Nova-ST",    9.27,        "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC002/bam_output",
  "Exp100-sampleSPC004",    "Mouse_liver",    "Nova-ST",    19.23,       "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp100-sampleSPC004/bam_output",
  "Exp65-sampleSPC014",     "Mouse_brain",    "Nova-ST",    20.27,       "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp65-sampleSPC014/bam_output",
  "Exp93-sampleSPC022",     "Human_lung",     "Nova-ST",    47.07,       "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/Exp93-sampleSPC022/bam_output",
  
  "spc002",     "Mouse_lung",     "Stereo-Seq",     13.47,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc002/bam_output",
  "spc004",     "Mouse_liver",    "Stereo-Seq",     17.21,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc004/bam_output",
  "spc014",     "Mouse_brain",    "Stereo-Seq",     22.09,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc014/bam_output",
  "spc022",     "Human_lung",     "Stereo-Seq",     44.24,   "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/spc022/bam_output",
  )

# Colors -------------------------------------------------------
tech_colors <- c(
  "Visium" = "#F4B18C",
  "Visium_HD_FFPE" = "#A6E18F",
  "Visium_HD_FF" = "#FFEC61",
  "BMKMANU_S1000" = "#D694F7",
  "BMKMANU_S3000" = "#F2A7E6",
  "Curio_Seeker" = "#9B9B9B",
  "Nova-ST" = "#F55151",
  "Stereo-Seq" = "#87DDFF"
)



sample_type_colors <- c(
  "Mouse_lung" = "#999999",
  "Mouse_liver" = "#E69F00",
  "Mouse_brain" = "#56B4E9",
  "Human_lung" = "#009E73",
  "Mouse_colon" = "#F0E442" # "#0072B2", "#D55E00", "#CC79A7"
)

# Output directories -----------------------------------------------------------------
tech_outdirs <- list(
  "Visium" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium/plots/bam_plots",
  "Visium_HD_FFPE" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/plots/FFPE/bam_plots",
  "Visium_HD_FF" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/plots/FF/bam_plots",
  "BMKMANU_S1000" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S1000/plots/bam_plots",
  "BMKMANU_S3000" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/plots/bam_plots",
  "Curio_Seeker" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/plots/bam_plots",
  "Nova-ST" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Nova-ST/plots/bam_plots",
  "Stereo-Seq" = "/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/Stereo-Seq/plots/bam_plots"
)

# Define function to read all metrics.csv for a given sample -------------------------
read_sample_metrics <- function(bam_dir, sample_name, sample_type, technology, area) {
  
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
      area = area
    )
}

# Load all samples --------------------------------------------------------
all_data <- map_dfr(1:nrow(sample_metadata), function(i) {
  row <- sample_metadata[i, ]
  read_sample_metrics(row$bam_dir, row$Sample, row$sample_type, row$Technology, row$area)
})
all_data <- all_data %>%
  mutate(
    distinct_per_mm2 = distinct / area,
    bound_pop_estimation_per_mm2 = bound_pop_estimation / area,
    bound_pop_lower_ci_per_mm2   = bound_pop_lower_ci / area,
    bound_pop_upper_ci_per_mm2   = bound_pop_upper_ci / area,
    sequencing_saturation_reads =
      if_else(corrected_reads > 0,
              (1 - (distinct / corrected_reads)) * 100,
              0),
    sequencing_saturation_UMIs =
      if_else(included_reads > 0,
              (1 - (distinct / included_reads)) * 100,
              0)
  )


# Function to make line plots ---------------------------------------------
make_line_plot <- function(df, xcol, ycol, outpath, color_by, palette, title, xlab, ylab) {
  df_plot <- df %>%
    rename(x_val = !!sym(xcol)) %>%
    rename(y_val = !!sym(ycol))
  
  use_area_corrected <- ycol == "distinct_per_mm2"
  
  percentage <- grepl(
    "percentage|sequencing_saturation", ycol)
  
  max_val <- max(df_plot$y_val, na.rm = TRUE)
  y_labels_func <- if (max_val > 1e6) {
    scales::label_number(scale = 1e-6, suffix = "M")
  } else {
    scales::comma
  }
  
  p <- ggplot(df_plot, aes(x = x_val, y = y_val,
                           color = !!sym(color_by), group = interaction(Sample, method))) +
    (if (!percentage) {
      ribbon_data <- df_plot %>%
        filter(method == "lc_extrap", x_val < 4e9) %>%
        mutate(
          ymin = if (use_area_corrected) lc_extrap_lower_ci / area else lc_extrap_lower_ci,
          ymax = if (use_area_corrected) lc_extrap_upper_ci / area else lc_extrap_upper_ci
        )
      
      # Smooth ymin and ymax using loess
      ribbon_data <- ribbon_data %>%
        group_by(Sample) %>%
        arrange(x_val) %>%
        mutate(
          ymin_smooth = predict(loess(ymin ~ x_val, span = 0.3)),
          ymax_smooth = predict(loess(ymax ~ x_val, span = 0.3))
        ) %>%
        ungroup()
      
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
      values = c(c_curve = "solid", lc_extrap = "dashed"),
      name = "Method"
    ) +
    (if (!percentage) {
      geom_vline(
        xintercept = 4e9,         # x-position of the line
        linetype = "dashed",      # dashed line
        color = "black",          # line color
        linewidth = 0.3                # line thickness
      ) 
    } else NULL) +
    (if (!percentage) {
      geom_point(
        data = df_plot %>%
          filter(method == "bound_pop") %>%
          mutate(x_val = 4.2e9, y_val = if (use_area_corrected) bound_pop_estimation_per_mm2 else bound_pop_estimation),
        aes(x = x_val, y = y_val, color = !!sym(color_by)),
        shape = 18,
        position = position_jitter(width = 0.075e9, height = 0),
        size = 2,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    } else NULL) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      text = element_text(color = "black", size = 12),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(linewidth=0.2, color="#DFDFDF"),
      panel.grid.major.y = element_line(linewidth=0.2, color="#DFDFDF"),
      axis.text = element_text(color = "black", size = 12),
      axis.title.x=element_text(size = 12, margin = margin(r = 5)),
      axis.title.y=element_blank(),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    labs(title = title, x = xlab, y = ylab, color = color_by) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) + 
    scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, suffix = "M"),
      breaks = pretty_breaks(n = 5),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(
      breaks = if (percentage) seq(0, 100, by = 10) else pretty_breaks(n = 6),
      labels = y_labels_func,
      limits = if (percentage) c(0, 100) else c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    guides(color = guide_legend(reverse = TRUE))
  
  ggsave(outpath,
         plot = p, width = 8, height = 4)
}

# Generate plots per sample -----------------------------------------------
samples <- distinct(all_data, Sample, bam_dir)
for (i in seq_len(nrow(samples))) {
  sub <- filter(all_data, Sample == samples$Sample[i])
  outdir <- file.path(samples$bam_dir[i], "plots", "preseq")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  sub_all <- filter(sub, scope == "all")
  make_line_plot(sub_all, "corrected_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_all_reads.png"), "Technology", tech_colors, "Number of Unique UMIs", "Reads", "Number of Unique UMIs")
  make_line_plot(sub_all, "included_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_UMI_reads.png"), "Technology", tech_colors, "Number of Unique UMIs", "Total UMIs", "Number of Unique UMIs")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_reads", file.path(outdir, "preseq_sequencing_saturation_reads_vs_all_reads.png"), "Technology", tech_colors, "Sequencing saturation - Reads", "Reads", "Sequencing saturation - Reads")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_all_reads.png"), "Technology", tech_colors, "Sequencing saturation - UMIs", "Reads", "Sequencing saturation - UMIs")
  make_line_plot(sub_all, "included_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_UMI_reads.png"), "Technology", tech_colors, "Sequencing saturation - UMIs", "Total UMIs", "Sequencing saturation - UMIs")
  
  sub_tissue <- filter(sub, scope == "tissue")
  make_line_plot(sub_tissue, "corrected_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_all_reads.png"), "Technology", tech_colors, "Number of Unique UMIs in tissue", "Reads", "Number of Unique UMIs in tissue")
  make_line_plot(sub_all, "included_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_UMI_reads.png"), "Technology", tech_colors, "Number of Unique UMIs in tissue", "Total UMIs", "Number of Unique UMIs in tissue")
  make_line_plot(sub_all, "corrected_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_all_reads.png"), "Technology", tech_colors, "Number of Unique UMIs in tissue per mm2", "Reads", "Number of Unique UMIs in tissue per mm2")
  make_line_plot(sub_all, "included_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_UMI_reads.png"), "Technology", tech_colors, "Number of Unique UMIs in tissue per mm2", "Total UMIs", "Number of Unique UMIs in tissue per mm2")
  }

# Generate plots per technology -------------------------------------------
tech_list <- unique(all_data$Technology)
for (tech in tech_list) {
  sub <- all_data[all_data$Technology == tech, ]
  
  outdir <- file.path(tech_outdirs[[tech]], "preseq")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  sub_all <- filter(sub, scope == "all")
  make_line_plot(sub_all, "corrected_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_all_reads.png"), "sample_type", sample_type_colors, "Number of Unique UMIs", "Reads", "Number of Unique UMIs")
  make_line_plot(sub_all, "included_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_UMI_reads.png"), "sample_type", sample_type_colors, "Number of Unique UMIs", "Total UMIs", "Number of Unique UMIs")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_reads", file.path(outdir, "preseq_sequencing_saturation_reads_vs_all_reads.png"), "sample_type", sample_type_colors, "Sequencing saturation - Reads", "Reads", "Sequencing saturation - Reads")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_all_reads.png"), "sample_type", sample_type_colors, "Sequencing saturation - UMIs", "Reads", "Sequencing saturation - UMIs")
  make_line_plot(sub_all, "included_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_UMI_reads.png"), "sample_type", sample_type_colors, "Sequencing saturation - UMIs", "Total UMIs", "Sequencing saturation - UMIs")
  
  sub_tissue <- filter(sub, scope == "tissue")
  make_line_plot(sub_tissue, "corrected_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_all_reads.png"), "sample_type", sample_type_colors, "Number of Unique UMIs in tissue", "Reads", "Number of Unique UMIs in tissue")
  make_line_plot(sub_all, "included_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_UMI_reads.png"), "sample_type", sample_type_colors, "Number of Unique UMIs in tissue", "Total UMIs", "Number of Unique UMIs in tissue")
  make_line_plot(sub_all, "corrected_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_all_reads.png"), "sample_type", sample_type_colors, "Number of Unique UMIs in tissue per mm2", "Reads", "Number of Unique UMIs in tissue per mm2")
  make_line_plot(sub_all, "included_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_UMI_reads.png"), "sample_type", sample_type_colors, "Number of Unique UMIs in tissue per mm2", "Total UMIs", "Number of Unique UMIs in tissue per mm2")
  
}

# Generate plots per sample type ------------------------------------------
sample_type_list <- unique(all_data$sample_type)
for (sample_type in sample_type_list) {
  sub <- all_data[all_data$sample_type == sample_type, ]
  
  outdir <- file.path("/data/groups/technologies/spatial.catalyst/Projects/2024-07-UTBenchmark-SpC/data/processed/plots/bam_plots", sample_type, "preseq")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  sub_all <- filter(sub, scope == "all")
  make_line_plot(sub_all, "corrected_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_all_reads.png"), "Technology", tech_colors, "Number of Unique UMIs", "Reads", "Number of Unique UMIs")
  make_line_plot(sub_all, "included_reads", "distinct", file.path(outdir, "preseq_all_UMIs_vs_UMI_reads.png"), "Technology", tech_colors, "Number of Unique UMIs", "Total UMIs", "Number of Unique UMIs")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_reads", file.path(outdir, "preseq_sequencing_saturation_reads_vs_all_reads.png"), "Technology", tech_colors, "Sequencing saturation - Reads", "Reads", "Sequencing saturation - Reads")
  make_line_plot(sub_all, "corrected_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_all_reads.png"), "Technology", tech_colors, "Sequencing saturation - UMIs", "Reads", "Sequencing saturation - UMIs")
  make_line_plot(sub_all, "included_reads", "sequencing_saturation_UMIs", file.path(outdir, "preseq_sequencing_saturation_UMIs_vs_UMI_reads.png"), "Technology", tech_colors, "Sequencing saturation - UMIs", "Total UMIs", "Sequencing saturation - UMIs")
  
  sub_tissue <- filter(sub, scope == "tissue")
  make_line_plot(sub_tissue, "corrected_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_all_reads.png"), "Technology", tech_colors, "Number of Unique UMIs in tissue", "Reads", "Number of Unique UMIs in tissue")
  make_line_plot(sub_all, "included_reads", "distinct", file.path(outdir, "preseq_UMIs_in_tissue_vs_UMI_reads.png"), "Technology", tech_colors, "Number of Unique UMIs in tissue", "Total UMIs", "Number of Unique UMIs in tissue")
  make_line_plot(sub_all, "corrected_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_all_reads.png"), "Technology", tech_colors, "Number of Unique UMIs in tissue per mm2", "Reads", "Number of Unique UMIs in tissue per mm2")
  make_line_plot(sub_all, "included_reads", "distinct_per_mm2", file.path(outdir, "preseq_UMIs_in_tissue_per_mm2_vs_UMI_reads.png"), "Technology", tech_colors, "Number of Unique UMIs in tissue per mm2", "Total UMIs", "Number of Unique UMIs in tissue per mm2")
  
  
}
