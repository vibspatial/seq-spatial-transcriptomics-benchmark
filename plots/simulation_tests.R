# ============================================================
# UMI Saturation Simulation
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# ---------------------------
# Function: simulate sampling from N unique objects
# ---------------------------
simulate_sampling <- function(N, max_reads = 10000, step = 100, seed = 123) {
  set.seed(seed)
  
  objects <- 1:N
  results <- list()
  sampled <- integer(0)
  
  for (R in seq(step, max_reads, by = step)) {
    draws <- sample(objects, size = step, replace = TRUE)
    sampled <- c(sampled, draws)
    
    unique_found <- length(unique(sampled))
    saturation <- (1 - (unique_found / R)) * 100
    
    results[[length(results) + 1]] <- tibble(
      total_reads = R,
      unique_found = unique_found,
      sequencing_saturation_reads = (1 - (unique_found / R)) * 100,
      N_library = N
    )
  }
  bind_rows(results)
}

# ============================================================
# 1. BASELINE: Varying library size
# ============================================================
all_sims <- bind_rows(lapply(seq(100, 1000, by = 100), function(N) {
  simulate_sampling(N, max_reads = 10000, step = 100)
}))

# ============================================================
# 2. Sampling efficiency loss (library size = 500)
# ============================================================
simulate_sampling_loss <- function(N, loss, max_reads = 10000, step = 100, seed = 123) {
  set.seed(seed)
  
  objects <- 1:N
  sampled <- integer(0)
  cumulative_usable_reads <- 0
  results <- list()
  
  for (R in seq(step, max_reads, by = step)) {
    # Draw 'step' reads
    draws <- sample(objects, size = step, replace = TRUE)
    
    # Apply loss: only a fraction of reads are usable
    usable <- sample(draws, size = round(step * (1 - loss)))
    
    # Track cumulative usable reads
    cumulative_usable_reads <- cumulative_usable_reads + length(usable)
    
    # Update sampled UMIs
    sampled <- c(sampled, usable)
    unique_found <- length(unique(sampled))
    
    # Store results
    results[[length(results) + 1]] <- tibble(
      total_reads = R,
      unique_found = unique_found,
      cumulative_usable_reads = cumulative_usable_reads,
      sequencing_saturation_reads = (1 - (unique_found / R)) * 100,
      sequencing_saturation_UMIs   = (1 - (unique_found / cumulative_usable_reads)) * 100,
      loss = paste0(loss * 100, "%")
    )
  }
  
  bind_rows(results)
}

# run across loss levels
base_500 <- simulate_sampling(500, max_reads = 10000, step = 100)
loss_levels <- seq(0, 0.9, by = 0.1)
loss_sims <- map_dfr(loss_levels, ~simulate_sampling_loss(500, .x))

# ============================================================
# 3. Scaling unique counts (library size = 500)
# ============================================================
scale_levels <- seq(0.1, 1, by = 0.1)
scaled_sims <- map_dfr(scale_levels, function(scale) {
  tibble(
    total_reads = base_500$total_reads,
    unique_found = base_500$unique_found * scale,
    sequencing_saturation_reads = (1 - ((base_500$unique_found * scale) / base_500$total_reads)) * 100,
    scale = paste0(scale * 100, "%")
  )
})

# ============================================================
# PLOTS
# ============================================================

# --- Baseline library size curves ---
p1 <- ggplot(all_sims, aes(x = total_reads, y = unique_found, color = factor(N_library))) +
  geom_line(size = 1.1) +
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
  labs(title = "Unique UMIs vs Reads (Different Library Sizes)",
       x = "Total Reads", y = "Unique UMIs Found", color = "Library Size (N)")

p2 <- ggplot(all_sims, aes(x = total_reads, y = sequencing_saturation_reads, color = factor(N_library))) +
  geom_line(size = 1.1) +
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
  labs(title = "Saturation vs Reads (Different Library Sizes)",
       x = "Total Reads", y = "Saturation - Reads (%)", color = "Library Size (N)")

# --- Sampling loss curves ---
p3 <- ggplot(loss_sims, aes(x = total_reads, y = unique_found, color = loss)) +
  geom_line(size = 1.1) +
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
  labs(title = "Unique UMIs vs Reads (Sampling Loss, N=500)",
       x = "Total Reads", y = "Unique UMIs Found", color = "Loss (%)")

p4 <- ggplot(loss_sims, aes(x = total_reads, y = sequencing_saturation_reads, color = loss)) +
  geom_line(size = 1.1) +
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
  labs(title = "Saturation vs Reads (Sampling Loss, N=500)",
       x = "Total Reads", y = "Saturation - Reads (%)", color = "Loss (%)")

p5 <- ggplot(loss_sims, aes(x = total_reads, y = sequencing_saturation_UMIs, color = loss)) +
  geom_line(size = 1.1) +
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
  labs(title = "Saturation vs Reads (Sampling Loss, N=500)",
       x = "Total Reads", y = "Saturation - UMIs (%)", color = "Loss (%)")

# --- Scaling unique counts curves ---
p6 <- ggplot(scaled_sims, aes(x = total_reads, y = unique_found, color = scale)) +
  geom_line(size = 1.1) +
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
  labs(title = "Unique UMIs vs Reads (Scaling Unique Counts, N=500)",
       x = "Total Reads", y = "Unique UMIs Found (scaled)", color = "Scaling (%)")

p7 <- ggplot(scaled_sims, aes(x = total_reads, y = sequencing_saturation_reads, color = scale)) +
  geom_line(size = 1.1) +
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
  labs(title = "Saturation vs Reads (Scaling Unique Counts, N=500)",
       x = "Total Reads", y = "Saturation - Reads (%)", color = "Scaling (%)")

# Save plots
output_dir <- "E:/2024-07-UTBenchmark-SpC/data/processed/plots/bam_plots/simulation_tests"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "p1_unique_vs_reads.png"), p1, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "p2_saturation_vs_reads.png"), p2, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "p3_loss_unique_vs_reads.png"), p3, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "p4_loss_saturation_reads.png"), p4, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "p5_loss_saturation_UMIs.png"), p5, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "p6_scaled_unique_vs_reads.png"), p6, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "p7_scaled_saturation_vs_reads.png"), p7, width = 7, height = 5, dpi = 300)

