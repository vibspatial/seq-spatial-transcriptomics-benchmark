# Working directory -------------------------------------------------------
setwd("E:/2024-07-UTBenchmark-SpC/data/processed/plots/correlation")

# Import libraries --------------------------------------------------------
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

# Read data ---------------------------------------------------------------
tissue <- 'human_lung'

df_Visium <- read.csv2("E:/2024-07-UTBenchmark-SpC/data/processed/Visium/SCA015/subsampled_100M/harpy/gene_total_counts.csv", header = TRUE, sep=',') %>%
  select(c(total_counts_corrected, gene)) %>%
  rename(Visium = total_counts_corrected)

df_VisiumHD_FFPE <- read.csv2("E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA010/subsampled_100M/harpy/gene_total_counts.csv", header = TRUE, sep=',') %>%
  select(c(total_counts_corrected, gene)) %>%
  rename(Visium_HD_FFPE = total_counts_corrected)

df_VisiumHD_FF <- read.csv2("E:/2024-07-UTBenchmark-SpC/data/processed/Visium_HD/SCA020/subsampled_100M/harpy/gene_total_counts.csv", header = TRUE, sep=',') %>%
  select(c(total_counts_corrected, gene)) %>%
  rename(Visium_HD_FF = total_counts_corrected)

df_S3000 <- read.csv2("E:/2024-07-UTBenchmark-SpC/data/processed/BMK_S3000/exp102/subsampled_100M/harpy/gene_total_counts.csv", header = TRUE, sep=',') %>%
  select(c(total_counts_corrected, gene)) %>%
  rename(BMKMANU_S3000 = total_counts_corrected)

df_Curioseeker <- read.csv2("E:/2024-07-UTBenchmark-SpC/data/processed/Curio_seeker/SC14/subsampled_100M/harpy/gene_total_counts.csv", header = TRUE, sep=',') %>%
  select(c(total_counts_corrected, gene)) %>%
  rename(Curio_Seeker = total_counts_corrected)

columns = c(
  "Visium",
  "Visium_HD_FFPE", 
  "Visium_HD_FF", 
  "BMKMANU_S3000", 
  "Curio_Seeker"
  )

# Combine genes -----------------------------------------------------------
genes <- rbind(select(df_Visium, gene),
               select(df_VisiumHD_FFPE, gene),
               select(df_VisiumHD_FF, gene),
               select(df_S3000, gene),
               select(df_Curioseeker, gene))
genes <- unique(genes)

# Append data -------------------------------------------------------------
correlation <- genes %>%
  left_join(df_Visium, by = "gene") %>%
  left_join(df_VisiumHD_FFPE, by = "gene") %>%
  left_join(df_VisiumHD_FF, by = "gene") %>%
  left_join(df_S3000, by = "gene") %>%
  left_join(df_Curioseeker, by = "gene") %>%
  mutate(across(columns, as.numeric)) %>%
  mutate(across(columns, ~ifelse(. <= 0, NA, .)))

# Raw correlation plot ---------------------------------------------------------
diagfun <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) + 
    geom_density(fill = "#5A2A81", alpha = 0.5, color = "black") + 
    scale_x_continuous(limits = c(1, 50000)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          text = element_text(color = "black", size = 6), 
          axis.text = element_text(size = 7, color = "black"), 
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
}

lowerfun <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "#AC94C0", size = 0.3, alpha = 0.6, na.rm = TRUE) +
    scale_x_continuous(limits = c(1, 50000)) +
    scale_y_continuous(limits = c(1, 50000)) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE, 
      color = "red",
      linetype = "dashed",
      linewidth = 0.6
    ) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth=0.6) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(color = "black", size = 6),
          axis.text = element_text(size = 7, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
}

upperfun <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) +
    stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", contour = TRUE, alpha = 0.6, na.rm = TRUE, n = 500, bins = 30) +
    scale_fill_viridis_c() +
    scale_x_continuous(limits = c(1, 50000)) +
    scale_y_continuous(limits = c(1, 50000)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth=0.6) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(color = "black", size = 6),
          axis.text = element_text(size = 7, color="black"),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
    stat_cor(
      aes(label = ..r.label..), 
      method = "pearson",
      label.x.npc = "left", label.y.npc = "top", size = 3, hjust = 0
    )
}

plot_correlation_raw <- ggpairs(
  data = correlation,
  columns = columns,
  diag = list(continuous = wrap(diagfun)),
  lower = list(continuous = wrap(lowerfun)),
  upper = list(continuous = wrap(upperfun))
) + theme(strip.text = element_text(size = 7))

ggsave(paste0(tissue, "_correlation_raw.png"),
       plot_correlation_raw, width = 7, height = 7, dpi = 300)

# Log correlation plot ---------------------------------------------------------
diagfun_log <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) + 
    geom_density(fill = "#5A2A81", alpha = 0.5, color = "black") + 
    scale_x_log10(limits = c(1, 50000)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          text = element_text(color = "black", size = 6),
          axis.text = element_text(size = 7, color = "black"), 
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
}

lowerfun_log <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "#AC94C0", size = 0.3, alpha = 0.6, na.rm = TRUE) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE, 
      color = "red",
      linetype = "dashed",
      linewidth = 0.6
    ) +
    scale_x_log10(limits = c(1, 50000)) +
    scale_y_log10(limits = c(1, 50000)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth=0.6) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(color = "black", size = 6),
          axis.text = element_text(size = 7, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
}

upperfun_log <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) +
    stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", contour = TRUE, alpha = 0.6, na.rm = TRUE, n = 500, bins = 30) +
    scale_fill_viridis_c() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth=0.6) +
    scale_x_log10(limits = c(1, 50000)) +
    scale_y_log10(limits = c(1, 50000)) +
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

plot_correlation_log <- ggpairs(
  data = correlation,
  columns = columns,
  diag = list(continuous = wrap(diagfun_log)),
  lower = list(continuous = wrap(lowerfun_log)),
  upper = list(continuous = wrap(upperfun_log))
) + theme(strip.text = element_text(size = 7))

ggsave(paste0(tissue, "_correlation_log.png"),
       plot_correlation_log, width = 7, height = 7, dpi = 300)