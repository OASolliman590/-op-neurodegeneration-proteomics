#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(yaml)
  library(scales)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, args)
  if (length(idx) == 0) return(getwd())
  normalizePath(dirname(sub(file_arg, "", args[idx][1])), mustWork = TRUE)
}

bool_present <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y", "t")
}

direction_from_fc <- function(fc, thr = 0.3) {
  ifelse(is.na(fc), "absent",
         ifelse(fc > thr, "up", ifelse(fc < -thr, "down", "flat")))
}

mk_label <- function(fc, marker = "") {
  ifelse(is.na(fc), "", sprintf("%.1f%s", fc, marker))
}

z_by_dataset <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(ifelse(is.na(x), NA_real_, 0))
  (x - m) / s
}

safe_max_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

main <- function() {
  root <- normalizePath(file.path(get_script_dir(), ".."), mustWork = FALSE)
  if (!file.exists(file.path(root, "config", "targets.yaml"))) {
    root <- normalizePath(getwd(), mustWork = TRUE)
  }
  out_fig <- file.path(root, "results", "figures")
  out_analysis <- file.path(root, "results", "analysis")
  dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_analysis, recursive = TRUE, showWarnings = FALSE)

  panel <- c("ACTG1", "DNAH9", "GPX3", "VWF", "C4B", "CD44", "CFHR2", "ITIH3", "LRG1", "MYH7B")

  # ------------------------------------------------------------------------
  # Load sources
  # ------------------------------------------------------------------------
  cfg <- yaml::read_yaml(file.path(root, "config", "targets.yaml"))
  targets <- cfg$targets %||% list()

  gene_full <- sapply(panel, function(g) {
    tg <- targets[[g]]
    if (is.null(tg)) g else (tg$full_name %||% g)
  }, USE.NAMES = FALSE)
  names(gene_full) <- panel
  gene_label_map <- setNames(paste0(panel, " - ", gene_full), panel)

  op_xlsx <- readxl::read_xlsx(file.path(root, "log fold change 10 Markers.xlsx"))
  op_df <- op_xlsx %>%
    transmute(
      gene = as.character(Marker),
      log2fc = as.numeric(`log2FC (Chronic/control)`)
    ) %>%
    filter(gene %in% panel) %>%
    mutate(
      dataset = "OP chronic",
      dataset_short = "OP chronic",
      section = "OP Exposure",
      disease = "OP",
      accession = "OP_DISCOVERY",
      n_case = NA_real_,
      n_control = NA_real_,
      present = !is.na(log2fc)
    ) %>%
    select(accession, disease, section, dataset, dataset_short, gene, present, n_case, n_control, log2fc)

  eprot <- readr::read_csv(file.path(root, "results", "eprot_query.csv"), show_col_types = FALSE)
  eprot_keep <- c("E-PROT-31", "E-PROT-61", "E-PROT-65")
  eprot_labels <- c(
    "E-PROT-31" = "AD brain (E-PROT-31)",
    "E-PROT-61" = "AD brain (E-PROT-61)",
    "E-PROT-65" = "PD brain (E-PROT-65)"
  )
  eprot_df <- eprot %>%
    filter(accession %in% eprot_keep, gene %in% panel) %>%
    mutate(
      section = "Brain Proteomics",
      dataset = eprot_labels[accession],
      dataset_short = dataset,
      n_case = suppressWarnings(as.numeric(n_case)),
      n_control = suppressWarnings(as.numeric(n_control)),
      present = bool_present(present),
      log2fc = suppressWarnings(as.numeric(log2fc))
    ) %>%
    select(accession, disease, section, dataset, dataset_short, gene, present, n_case, n_control, log2fc)

  pride_old <- readr::read_csv(file.path(root, "results", "pride_query.csv"), show_col_types = FALSE)
  ad_csf_df <- pride_old %>%
    filter(accession == "PXD016278", gene %in% panel) %>%
    mutate(
      section = "CSF Proteomics",
      dataset = "AD CSF (PXD016278)",
      dataset_short = dataset,
      n_case = suppressWarnings(as.numeric(n_case)),
      n_control = suppressWarnings(as.numeric(n_control)),
      mean_case = suppressWarnings(as.numeric(mean_case)),
      mean_ctrl = suppressWarnings(as.numeric(mean_ctrl)),
      log2fc = ifelse(!is.na(mean_case) & !is.na(mean_ctrl), mean_case - mean_ctrl, NA_real_),
      present = bool_present(present)
    ) %>%
    select(accession, disease, section, dataset, dataset_short, gene, present, n_case, n_control, log2fc)

  pride_new <- readr::read_csv(file.path(root, "results", "pride_discovery", "pride_ad_pd_ms_query.csv"), show_col_types = FALSE)
  pride_keep <- c("PXD026491", "PXD045058")
  pride_labels <- c(
    "PXD026491" = "PD CSF (PXD026491)",
    "PXD045058" = "MS CSF (PXD045058)"
  )
  pride_df <- pride_new %>%
    filter(accession %in% pride_keep, gene %in% panel) %>%
    mutate(
      section = "CSF Proteomics",
      dataset = pride_labels[accession],
      dataset_short = dataset,
      n_case = suppressWarnings(as.numeric(n_case)),
      n_control = suppressWarnings(as.numeric(n_control)),
      log2fc = suppressWarnings(as.numeric(log2fc)),
      present = bool_present(present)
    ) %>%
    select(accession, disease, section, dataset, dataset_short, gene, present, n_case, n_control, log2fc)

  combined <- bind_rows(op_df, eprot_df, ad_csf_df, pride_df)

  dataset_order <- c(
    "OP chronic",
    "AD brain (E-PROT-31)",
    "AD brain (E-PROT-61)",
    "PD brain (E-PROT-65)",
    "AD CSF (PXD016278)",
    "PD CSF (PXD026491)",
    "MS CSF (PXD045058)"
  )
  section_order <- c("OP Exposure", "Brain Proteomics", "CSF Proteomics")

  # Ensure complete panel x dataset matrix
  combined <- combined %>%
    mutate(
      dataset = factor(dataset, levels = dataset_order),
      section = factor(section, levels = section_order),
      gene = as.character(gene)
    ) %>%
    tidyr::complete(
      dataset = factor(dataset_order, levels = dataset_order),
      gene = panel,
      fill = list(
        accession = NA_character_,
        disease = NA_character_,
        section = NA_character_,
        dataset_short = NA_character_,
        present = FALSE,
        n_case = NA_real_,
        n_control = NA_real_,
        log2fc = NA_real_
      )
    ) %>%
    group_by(dataset) %>%
    fill(section, .direction = "downup") %>%
    ungroup() %>%
    mutate(
      section = factor(section, levels = section_order),
      disease = coalesce(disease, "NA")
    )

  # Dataset metadata for labels
  meta <- combined %>%
    group_by(dataset, section) %>%
    summarise(
      n_case = safe_max_num(n_case),
      n_control = safe_max_num(n_control),
      .groups = "drop"
    ) %>%
    mutate(
      dataset_label = ifelse(
        !is.na(n_case) & !is.na(n_control) & n_case > 0 & n_control > 0,
        sprintf("%s\n(n=%s/%s)", as.character(dataset), format(n_case, trim = TRUE), format(n_control, trim = TRUE)),
        as.character(dataset)
      )
    )

  combined <- combined %>%
    left_join(meta %>% select(dataset, dataset_label), by = "dataset")

  # OP direction reference (from experimental OP table)
  op_ref <- combined %>%
    filter(section == "OP Exposure") %>%
    transmute(gene, op_direction = direction_from_fc(log2fc))

  combined <- combined %>%
    left_join(op_ref, by = "gene") %>%
    mutate(
      direction = direction_from_fc(log2fc),
      concordance = case_when(
        section == "OP Exposure" ~ "ref",
        direction %in% c("up", "down") & op_direction %in% c("up", "down") & direction == op_direction ~ "yes",
        direction %in% c("up", "down") & op_direction %in% c("up", "down") & direction != op_direction ~ "no",
        TRUE ~ "na"
      ),
      marker = case_when(
        concordance == "yes" ~ " [C]",
        concordance == "no" ~ " [D]",
        TRUE ~ ""
      ),
      label = mk_label(log2fc, marker),
      zscore = ifelse(present, log2fc, NA_real_)
    ) %>%
    group_by(dataset) %>%
    mutate(zscore = z_by_dataset(zscore)) %>%
    ungroup() %>%
    mutate(
      gene_label = factor(gene_label_map[gene], levels = rev(gene_label_map[panel])),
      dataset_label = factor(dataset_label, levels = meta$dataset_label[match(dataset_order, meta$dataset)]),
      section = factor(section, levels = section_order)
    )

  # ------------------------------------------------------------------------
  # Plot A: Effect-size heatmap (dataset-standardized color, raw log2FC labels)
  # ------------------------------------------------------------------------
  p_effect <- ggplot(combined, aes(x = dataset_label, y = gene_label, fill = zscore)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 3.0, color = "black", na.rm = TRUE) +
    facet_grid(. ~ section, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, na.value = "#e9e9e9",
      limits = c(-2.5, 2.5), oob = scales::squish,
      name = "Within-dataset\nz-score"
    ) +
    labs(
      title = "OP Signature vs Neurodegeneration Proteomics",
      subtitle = "Tile color: within-dataset standardized effect. Tile text: raw log2FC; [C]=concordant, [D]=discordant to OP direction.",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    )

  # ------------------------------------------------------------------------
  # Plot B: Concordance map (explicit validity view)
  # ------------------------------------------------------------------------
  concordance_levels <- c("ref", "yes", "no", "na")
  combined$concordance <- factor(combined$concordance, levels = concordance_levels)

  p_conc <- ggplot(combined, aes(x = dataset_label, y = gene_label, fill = concordance)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(log2fc), "", sprintf("%.1f", log2fc))), size = 2.8, color = "black", na.rm = TRUE) +
    facet_grid(. ~ section, scales = "free_x", space = "free_x") +
    scale_fill_manual(
      values = c(
        ref = "#1f3b73",
        yes = "#2e7d32",
        no = "#c62828",
        na = "#d9d9d9"
      ),
      drop = FALSE,
      name = "Concordance"
    ) +
    labs(
      title = "Concordance to OP Direction",
      subtitle = "Green=yes, Red=no, Gray=not scoreable (flat/absent). Values are raw log2FC.",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    )

  final_plot <- p_effect / p_conc + plot_layout(heights = c(1.08, 0.92))

  fig_path_png <- file.path(out_fig, "grand_finale_op_neurodegeneration_ggplot.png")
  fig_path_pdf <- file.path(out_fig, "grand_finale_op_neurodegeneration_ggplot.pdf")
  ggsave(fig_path_png, final_plot, width = 17, height = 12, dpi = 320, bg = "white")
  ggsave(fig_path_pdf, final_plot, width = 17, height = 12, bg = "white")

  # ------------------------------------------------------------------------
  # Scientific validity appraisal (data-driven)
  # ------------------------------------------------------------------------
  non_op <- combined %>% filter(section != "OP Exposure")
  scoreable <- non_op %>% filter(concordance %in% c("yes", "no"))
  yes_n <- sum(scoreable$concordance == "yes", na.rm = TRUE)
  no_n <- sum(scoreable$concordance == "no", na.rm = TRUE)
  na_n <- sum(non_op$concordance == "na", na.rm = TRUE)

  ds_validity <- non_op %>%
    group_by(dataset, section) %>%
    summarise(
      n_case = safe_max_num(n_case),
      n_control = safe_max_num(n_control),
      n_present = sum(present & !is.na(log2fc), na.rm = TRUE),
      n_scoreable = sum(concordance %in% c("yes", "no"), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      comparability_tier = case_when(
        !is.na(n_case) & !is.na(n_control) & n_case > 0 & n_control > 0 ~ "Tier A (usable case-control metadata)",
        TRUE ~ "Tier C (insufficient case-control metadata)"
      )
    )

  cfg_dir <- sapply(panel, function(g) tolower((targets[[g]]$chronic_direction %||% "none")), USE.NAMES = TRUE)
  op_dir <- combined %>%
    filter(section == "OP Exposure") %>%
    transmute(gene, op_dir = direction_from_fc(log2fc))
  op_dir_vec <- setNames(op_dir$op_dir, op_dir$gene)
  mismatch <- panel[cfg_dir[panel] %in% c("up", "down") & op_dir_vec[panel] %in% c("up", "down") & cfg_dir[panel] != op_dir_vec[panel]]

  data_export <- combined %>%
    select(accession, disease, section, dataset, dataset_label, gene, gene_label, present, n_case, n_control, log2fc, direction, op_direction, concordance, zscore)
  write_csv(data_export, file.path(out_analysis, "grand_finale_ggplot_data.csv"))
  write_csv(ds_validity, file.path(out_analysis, "grand_finale_dataset_validity.csv"))

  appraisal_path <- file.path(out_analysis, "scientific_validity_appraisal.md")
  lines <- c(
    "# Scientific Validity Appraisal (Figure-Level)",
    "",
    sprintf("Date: %s", as.character(Sys.Date())),
    "",
    "## 1) What Is Scientifically Strong",
    sprintf("- All comparisons are within proteomics modality, reducing cross-platform noise."),
    sprintf("- Figure includes %d non-OP datasets and %d total scoreable gene-dataset cells.", n_distinct(non_op$dataset), nrow(scoreable)),
    sprintf("- Concordance among scoreable cells: **yes=%d**, **no=%d**.", yes_n, no_n),
    "",
    "## 2) Critical Validity Risks",
    sprintf("- **Direction-definition conflict**: OP sign from `log fold change 10 Markers.xlsx` disagrees with `targets.yaml` chronic directions for: %s.",
            ifelse(length(mismatch) > 0, paste(mismatch, collapse = ", "), "none")),
    "- Brain and CSF are biologically distinct compartments; agreement should be interpreted as cross-compartment convergence, not direct replication.",
    "- Some datasets still lack explicit case/control counts, limiting weighted inference.",
    "- Current visualization uses point estimates (log2FC) without confidence intervals or adjusted p-values.",
    "",
    "## 3) Dataset Comparability Tiers",
    "",
    "| Dataset | Section | n_case | n_control | n_present | n_scoreable | Tier |",
    "|---|---|---:|---:|---:|---:|---|"
  )
  for (i in seq_len(nrow(ds_validity))) {
    r <- ds_validity[i, ]
    lines <- c(
      lines,
      sprintf("| %s | %s | %s | %s | %d | %d | %s |",
              as.character(r$dataset),
              as.character(r$section),
              ifelse(is.na(r$n_case), "NA", as.character(as.integer(r$n_case))),
              ifelse(is.na(r$n_control), "NA", as.character(as.integer(r$n_control))),
              as.integer(r$n_present),
              as.integer(r$n_scoreable),
              as.character(r$comparability_tier))
    )
  }
  lines <- c(
    lines,
    "",
    "## 4) Methodology Guardrails for Manuscript Claims",
    "- Use spreadsheet-derived OP direction as the source-of-truth (`results/analysis/op_direction_source_of_truth.csv`).",
    "- Base primary inference on stratified Tier A analyses (brain and CSF reported separately).",
    "- Report Tier A + Tier B sensitivity analyses while excluding Tier C from primary inference.",
    "- Frame findings as cross-disease overlap / partial concordance unless validated in independent OP-exposed human cohorts.",
    "",
    "## 5) Outputs Produced",
    sprintf("- Figure PNG: `%s`", fig_path_png),
    sprintf("- Figure PDF: `%s`", fig_path_pdf),
    sprintf("- Appraisal table: `%s`", file.path(out_analysis, "grand_finale_dataset_validity.csv")),
    sprintf("- Figure data table: `%s`", file.path(out_analysis, "grand_finale_ggplot_data.csv"))
  )
  writeLines(lines, appraisal_path)

  message("Saved figure: ", fig_path_png)
  message("Saved figure: ", fig_path_pdf)
  message("Saved appraisal: ", appraisal_path)
}

main()
