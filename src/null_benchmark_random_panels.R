#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  flag <- "--file="
  idx <- grep(flag, args)
  if (length(idx) == 0) return(getwd())
  normalizePath(dirname(sub(flag, "", args[idx][1])), mustWork = TRUE)
}

direction_from_fc <- function(fc, thr = 0.3) {
  ifelse(is.na(fc), "no_data",
         ifelse(fc > thr, "up", ifelse(fc < -thr, "down", "flat")))
}

concordance_rate <- function(direction, expected) {
  scoreable <- direction %in% c("up", "down") & expected %in% c("up", "down")
  n <- sum(scoreable, na.rm = TRUE)
  if (n == 0) return(list(rate = NA_real_, n_scoreable = 0L, n_yes = 0L, n_no = 0L))
  yes <- sum(direction[scoreable] == expected[scoreable], na.rm = TRUE)
  no <- n - yes
  list(rate = yes / n, n_scoreable = as.integer(n), n_yes = as.integer(yes), n_no = as.integer(no))
}

main <- function() {
  root <- normalizePath(file.path(get_script_dir(), ".."), mustWork = FALSE)
  if (!file.exists(file.path(root, "results", "analysis", "publication_effects_with_uncertainty.csv"))) {
    root <- normalizePath(getwd(), mustWork = TRUE)
  }

  out_analysis <- file.path(root, "results", "analysis")
  out_fig <- file.path(root, "results", "figures")
  dir.create(out_analysis, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

  canonical_path <- file.path(out_analysis, "publication_effects_with_uncertainty.csv")
  op_path <- file.path(out_analysis, "op_direction_source_of_truth.csv")
  matrix_path <- file.path(root, "results", "pride_discovery", "PXD045058_unique_genes_matrix.tsv")
  anno_path <- file.path(root, "results", "pride_discovery", "PXD045058_sample_annotations.xlsx")

  if (!file.exists(canonical_path) || !file.exists(op_path)) {
    stop("Missing canonical publication tables. Run publication_overlap_analysis.R first.")
  }
  if (!file.exists(matrix_path) || !file.exists(anno_path)) {
    stop("Missing PXD045058 full matrix files for null benchmark.")
  }

  canonical <- read_csv(canonical_path, show_col_types = FALSE)
  op <- read_csv(op_path, show_col_types = FALSE) %>% select(gene, op_direction)

  panel <- op$gene
  expected_panel <- setNames(op$op_direction, op$gene)

  obs <- canonical %>%
    filter(accession == "PXD045058", gene %in% panel) %>%
    distinct(gene, .keep_all = TRUE) %>%
    mutate(
      direction = case_when(
        effect_status %in% c("up", "down", "flat") ~ effect_status,
        TRUE ~ direction_from_fc(log2fc, thr = 0.3)
      ),
      expected = expected_panel[gene]
    )

  obs_metrics <- concordance_rate(obs$direction, obs$expected)

  anno <- read_excel(anno_path)
  if (!all(c("Run", "Diagnosis_group") %in% names(anno))) {
    stop("PXD045058 annotations missing required columns: Run, Diagnosis_group")
  }

  mat <- read.delim(matrix_path, check.names = FALSE)
  if (ncol(mat) < 3) stop("PXD045058 matrix has insufficient columns")
  colnames(mat)[1] <- "gene"

  sample_cols <- colnames(mat)[-1]
  run_key <- basename(sample_cols)
  run_key <- sub("\\.dia$", "", run_key)

  anno <- anno %>% mutate(Run = as.character(Run), Diagnosis_group = as.character(Diagnosis_group))
  case_runs <- anno %>% filter(Diagnosis_group == "MS") %>% pull(Run)
  ctrl_runs <- anno %>%
    filter(!(Diagnosis_group %in% c("MS", "QCpool_PILOT1", "QCpool_PILOT2"))) %>%
    pull(Run)

  case_idx <- which(run_key %in% case_runs)
  ctrl_idx <- which(run_key %in% ctrl_runs)
  if (length(case_idx) < 10 || length(ctrl_idx) < 10) {
    stop(sprintf("Insufficient mapped samples in matrix. case=%d control=%d", length(case_idx), length(ctrl_idx)))
  }

  x <- as.matrix(mat[, -1, drop = FALSE])
  suppressWarnings(mode(x) <- "numeric")
  mean_case <- rowMeans(x[, case_idx, drop = FALSE], na.rm = TRUE)
  mean_ctrl <- rowMeans(x[, ctrl_idx, drop = FALSE], na.rm = TRUE)

  universe <- tibble(
    gene = as.character(mat$gene),
    mean_case = mean_case,
    mean_ctrl = mean_ctrl,
    log2fc = mean_case - mean_ctrl,
    direction = direction_from_fc(log2fc, thr = 0.3)
  ) %>%
    filter(!is.na(gene), gene != "") %>%
    distinct(gene, .keep_all = TRUE)

  universe_det <- universe %>% filter(direction %in% c("up", "down", "flat"))
  if (nrow(universe_det) < length(panel) + 50) {
    stop("Gene universe too small for stable random-panel benchmark")
  }

  set.seed(47)
  n_sim <- 10000L
  expected_vec <- unname(expected_panel[panel])
  expected_vec <- expected_vec[expected_vec %in% c("up", "down")]
  if (length(expected_vec) == 0) stop("No up/down expected directions available in OP panel")

  sim_rows <- vector("list", n_sim)
  for (i in seq_len(n_sim)) {
    sel <- sample(universe_det$gene, size = length(panel), replace = FALSE)
    obs_dir <- universe_det$direction[match(sel, universe_det$gene)]
    exp_dir <- sample(expected_vec, size = length(panel), replace = TRUE)

    met <- concordance_rate(obs_dir, exp_dir)
    sim_rows[[i]] <- tibble(
      sim_id = i,
      concordance_rate = met$rate,
      n_scoreable = met$n_scoreable,
      n_yes = met$n_yes,
      n_no = met$n_no
    )
  }

  sim_tbl <- bind_rows(sim_rows)
  sim_tbl <- sim_tbl %>% mutate(matched_detectability = n_scoreable == obs_metrics$n_scoreable)

  matched_tbl <- sim_tbl %>% filter(matched_detectability)
  base_tbl <- if (nrow(matched_tbl) >= 200) matched_tbl else sim_tbl
  benchmark_mode <- if (nrow(matched_tbl) >= 200) "matched_detectability" else "all_random_panels"

  obs_rate <- obs_metrics$rate
  p_empirical <- mean(base_tbl$concordance_rate >= obs_rate, na.rm = TRUE)
  pct <- mean(base_tbl$concordance_rate <= obs_rate, na.rm = TRUE)

  summary_tbl <- tibble(
    dataset = "PXD045058",
    benchmark_mode = benchmark_mode,
    n_sim_total = nrow(sim_tbl),
    n_sim_used = nrow(base_tbl),
    observed_concordance_rate = obs_rate,
    observed_n_scoreable = obs_metrics$n_scoreable,
    observed_n_yes = obs_metrics$n_yes,
    observed_n_no = obs_metrics$n_no,
    null_mean_rate = mean(base_tbl$concordance_rate, na.rm = TRUE),
    null_sd_rate = sd(base_tbl$concordance_rate, na.rm = TRUE),
    empirical_p_ge_observed = p_empirical,
    empirical_percentile = pct,
    n_case = length(case_idx),
    n_control = length(ctrl_idx)
  )

  panel_obs_tbl <- obs %>%
    select(gene, log2fc, effect_status, op_direction) %>%
    mutate(dataset = "PXD045058")

  write_csv(sim_tbl, file.path(out_analysis, "null_benchmark_random_panels.csv"))
  write_csv(summary_tbl, file.path(out_analysis, "null_benchmark_summary.csv"))
  write_csv(panel_obs_tbl, file.path(out_analysis, "null_benchmark_panel_observed.csv"))

  p <- ggplot(base_tbl, aes(x = concordance_rate)) +
    geom_histogram(binwidth = 0.05, fill = "#4c78a8", color = "white", alpha = 0.9) +
    geom_vline(xintercept = obs_rate, color = "#e45756", linewidth = 1.1) +
    annotate("text", x = obs_rate, y = Inf, label = sprintf("Observed = %.2f", obs_rate),
             color = "#e45756", vjust = 1.8, hjust = -0.05, size = 3.3) +
    labs(
      title = "Random-panel null benchmark (PXD045058)",
      subtitle = sprintf("Mode: %s | empirical p = %.4f", benchmark_mode, p_empirical),
      x = "Concordance rate",
      y = "Random panel count"
    ) +
    theme_minimal(base_size = 11)

  ggsave(file.path(out_fig, "null_benchmark_random_panels.png"), p,
         width = 9, height = 5.6, dpi = 300, bg = "white")

  message("Saved null benchmark outputs:")
  message(" - ", file.path(out_analysis, "null_benchmark_summary.csv"))
  message(" - ", file.path(out_analysis, "null_benchmark_random_panels.csv"))
  message(" - ", file.path(out_fig, "null_benchmark_random_panels.png"))
}

main()
