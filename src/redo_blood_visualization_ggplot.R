#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(yaml)
  library(patchwork)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  flag <- "--file="
  idx <- grep(flag, args)
  if (length(idx) == 0) return(getwd())
  normalizePath(dirname(sub(flag, "", args[idx][1])), mustWork = TRUE)
}

main <- function() {
  root <- normalizePath(file.path(get_script_dir(), ".."), mustWork = FALSE)
  if (!file.exists(file.path(root, "config", "targets.yaml"))) {
    root <- normalizePath(getwd(), mustWork = TRUE)
  }

  t7_root <- "/Volumes/T7/5-Alzhimers_Parkisons_MS_External_Valida/op_external_validation_data"
  canonical_res <- file.path(t7_root, "results", "blood_validation_live_primary_pxd")
  local_res <- file.path(root, "results", "blood_validation_live_primary_pxd")
  input_res <- if (dir.exists(canonical_res)) canonical_res else local_res

  out_fig <- file.path(root, "results", "figures")
  out_analysis <- file.path(root, "results", "analysis")
  dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_analysis, recursive = TRUE, showWarnings = FALSE)

  fx_path <- file.path(input_res, "per_dataset_effects.csv")
  qc_path <- file.path(input_res, "per_dataset_qc.csv")
  if (!file.exists(fx_path) || !file.exists(qc_path)) {
    stop("Missing per_dataset_effects.csv or per_dataset_qc.csv in: ", input_res)
  }

  fx <- read_csv(fx_path, show_col_types = FALSE)
  qc <- read_csv(qc_path, show_col_types = FALSE)

  cfg <- yaml::read_yaml(file.path(root, "config", "targets.yaml"))
  panel <- names(cfg$targets)
  gene_full <- sapply(panel, function(g) cfg$targets[[g]]$full_name %||% g, USE.NAMES = FALSE)
  names(gene_full) <- panel
  gene_levels <- rev(paste0(panel, " - ", gene_full[panel]))

  disease_levels <- c("Alzheimers", "Parkinsons", "MS")
  status_levels <- c("ok", "blocked_no_group_labels", "blocked_no_quant_table", "blocked_no_id_mapping", "parse_failed")

  # Dataset metadata
  meta <- qc %>%
    mutate(
      disease = factor(disease, levels = disease_levels),
      status = factor(status, levels = status_levels),
      dataset_label = sprintf("%s\n%s | n=%s/%s",
                              accession, biospecimen,
                              ifelse(is.na(n_case), "0", as.character(n_case)),
                              ifelse(is.na(n_control), "0", as.character(n_control)))
    ) %>%
    arrange(disease, accession) %>%
    mutate(dataset_label = factor(dataset_label, levels = dataset_label))

  ok_labels <- meta %>% filter(status == "ok") %>% pull(dataset_label) %>% as.character()

  fig_df <- fx %>%
    left_join(meta %>% select(accession, disease, biospecimen, status, dataset_label), by = c("accession", "disease", "biospecimen")) %>%
    mutate(
      gene_label = factor(paste0(gene, " - ", gene_full[gene]), levels = gene_levels),
      present = tolower(as.character(present)) %in% c("true", "1", "yes", "y", "t"),
      log2fc_num = suppressWarnings(as.numeric(log2fc)),
      concordance = case_when(
        tolower(as.character(concordant)) == "yes" ~ "yes",
        tolower(as.character(concordant)) == "no" ~ "no",
        TRUE ~ "na"
      ),
      dir_pair = ifelse(is.na(direction) | is.na(expected), "", paste0(direction, " / ", expected))
    ) %>%
    filter(gene %in% panel)

  fig_ok <- fig_df %>% filter(as.character(dataset_label) %in% ok_labels)

  # Panel A: effect-size heatmap
  p_effect <- ggplot(fig_ok, aes(x = dataset_label, y = gene_label, fill = log2fc_num)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(present & !is.na(log2fc_num), sprintf("%.2f", log2fc_num), "")),
              size = 2.9, color = "black") +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, na.value = "#d9d9d9", limits = c(-3, 3), oob = scales::squish,
      name = "log2FC"
    ) +
    labs(
      title = "A) Blood Proteomics Effect Size (Analyzable Datasets Only)",
      subtitle = "Tiles show case-control log2FC per OP-panel gene. Grey = not detected / unavailable.",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "right"
    )

  # Panel B: concordance heatmap
  p_conc <- ggplot(fig_ok, aes(x = dataset_label, y = gene_label, fill = concordance)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(present, dir_pair, "")), size = 2.4, color = "black") +
    scale_fill_manual(
      values = c("yes" = "#2e7d32", "no" = "#c62828", "na" = "#d9d9d9"),
      breaks = c("yes", "no", "na"),
      name = "Concordance\nvs OP"
    ) +
    labs(
      title = "B) Directional Concordance",
      subtitle = "Cell text = observed / expected direction.",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "right"
    )

  # Panel C: QC / blockers map for all screened datasets
  qc_plot_df <- meta %>%
    mutate(
      dataset_qc_label = sprintf("%s | %s | n=%s/%s",
                                 accession, biospecimen,
                                 ifelse(is.na(n_case), "0", as.character(n_case)),
                                 ifelse(is.na(n_control), "0", as.character(n_control))),
      dataset_qc_label = factor(dataset_qc_label, levels = rev(dataset_qc_label)),
      qc_x = "dataset_status"
    )

  p_qc <- ggplot(qc_plot_df, aes(x = qc_x, y = dataset_qc_label, fill = status)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("genes=%s", n_panel_genes_present)), size = 2.4, color = "black") +
    scale_fill_manual(
      values = c(
        "ok" = "#2e7d32",
        "blocked_no_group_labels" = "#f9a825",
        "blocked_no_quant_table" = "#8e8e8e",
        "blocked_no_id_mapping" = "#6d4c41",
        "parse_failed" = "#4e342e"
      ),
      drop = FALSE,
      name = "QC status"
    ) +
    facet_grid(disease ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "C) Dataset QC and Blockers (All 21 Screened Datasets)",
      subtitle = "Tile text = number of OP-panel genes present in parsed table.",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.text.y = element_text(face = "bold", angle = 0),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "right"
    )

  final_plot <- (p_effect | p_conc) / p_qc +
    plot_layout(heights = c(1.0, 1.05), widths = c(1, 1)) +
    plot_annotation(
      title = "External Blood Validation: OP 10-Gene Panel (AD / PD / MS)",
      subtitle = sprintf("Input: %s | Analyzable datasets: %s/%s",
                         basename(input_res),
                         sum(qc$status == "ok", na.rm = TRUE),
                         nrow(qc)),
      theme = theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 10)
      )
    )

  png_out <- file.path(out_fig, "blood_validation_live_primary_ggplot.png")
  pdf_out <- file.path(out_fig, "blood_validation_live_primary_ggplot.pdf")
  ggsave(png_out, final_plot, width = 20, height = 14, dpi = 300, bg = "white")
  ggsave(pdf_out, final_plot, width = 20, height = 14, bg = "white")

  # Save figure-ready data for transparency
  write_csv(fig_ok, file.path(out_analysis, "blood_validation_live_primary_figure_data.csv"))
  write_csv(qc_plot_df, file.path(out_analysis, "blood_validation_live_primary_qc_plot_data.csv"))

  # Optional copy to T7 figure folder
  t7_fig <- file.path(t7_root, "results", "figures")
  if (dir.exists(file.path(t7_root, "results"))) {
    dir.create(t7_fig, recursive = TRUE, showWarnings = FALSE)
    file.copy(png_out, file.path(t7_fig, basename(png_out)), overwrite = TRUE)
    file.copy(pdf_out, file.path(t7_fig, basename(pdf_out)), overwrite = TRUE)
  }

  message("Saved: ", png_out)
  message("Saved: ", pdf_out)
}

main()

