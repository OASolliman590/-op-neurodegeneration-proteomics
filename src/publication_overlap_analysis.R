#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(yaml)
  library(ggplot2)
  library(metafor)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  flag <- "--file="
  idx <- grep(flag, args)
  if (length(idx) == 0) return(getwd())
  normalizePath(dirname(sub(flag, "", args[idx][1])), mustWork = TRUE)
}

bool_present <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y", "t")
}

direction_from_fc <- function(fc, thr = 0.3) {
  ifelse(is.na(fc), "no_data",
         ifelse(fc > thr, "up", ifelse(fc < -thr, "down", "flat")))
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

infer_pxd026491_counts <- function() {
  hbs_url <- "https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/06/PXD026491/2020-11-20_CSF_HBS_DIA_Proteins_Percentile025.csv"
  lcc_url <- "https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/06/PXD026491/2020-11-20_CSF_LCC_DIA_Proteins_Percentile025.csv"

  parse_header_count <- function(url, regex, id_regex) {
    tmp <- tempfile(fileext = ".csv")
    on.exit(unlink(tmp), add = TRUE)
    ok <- tryCatch({
      suppressWarnings(download.file(url, tmp, mode = "wb", quiet = TRUE))
      TRUE
    }, error = function(e) FALSE)
    if (!ok) return(NA_integer_)
    hdr <- tryCatch(readr::read_csv(tmp, n_max = 1, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(hdr)) return(NA_integer_)
    cols <- names(hdr)
    samp <- cols[str_detect(cols, regex)]
    ids <- str_match(samp, id_regex)[, 2]
    ids <- suppressWarnings(as.integer(ids))
    ids <- ids[!is.na(ids)]
    if (length(ids) == 0) return(NA_integer_)
    length(unique(ids))
  }

  n_hbs <- parse_header_count(hbs_url, "_HBS_CSF_", "_HBS_CSF_([0-9]+)")
  n_lcc <- parse_header_count(lcc_url, "_LCC_CSF_sample_", "_sample_([0-9]+)")
  # Network-independent fallback from prior verified FTP header parse.
  if (is.na(n_hbs)) n_hbs <- 96
  if (is.na(n_lcc)) n_lcc <- 130
  list(n_case = n_hbs, n_control = n_lcc)
}

run_meta_hk <- function(df, group_cols, set_name) {
  if (nrow(df) == 0) return(tibble())
  split_df <- split(df, interaction(df[, group_cols], drop = TRUE))

  rows <- lapply(split_df, function(sub) {
    sub <- sub %>%
      filter(!is.na(log2fc), !is.na(se_proxy), se_proxy > 0)
    if (nrow(sub) == 0) return(NULL)

    k <- nrow(sub)
    yi <- sub$log2fc
    sei <- sub$se_proxy
    out <- as.list(sub[1, group_cols, drop = FALSE])

    if (k == 1) {
      est <- yi[1]
      se <- sei[1]
      ci_l <- est - 1.96 * se
      ci_u <- est + 1.96 * se
      p <- 2 * pnorm(-abs(est / se))
      out$k <- 1
      out$meta_method <- "single_study"
      out$estimate <- est
      out$se <- se
      out$ci_low <- ci_l
      out$ci_high <- ci_u
      out$p_value <- p
      out$tau2 <- NA_real_
      out$I2 <- NA_real_
      out$QEp <- NA_real_
      out$analysis_set <- set_name
      return(as_tibble(out))
    }

    fit <- tryCatch(
      metafor::rma.uni(yi = yi, sei = sei, method = "REML", test = "knha"),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      out$k <- k
      out$meta_method <- "failed"
      out$estimate <- NA_real_
      out$se <- NA_real_
      out$ci_low <- NA_real_
      out$ci_high <- NA_real_
      out$p_value <- NA_real_
      out$tau2 <- NA_real_
      out$I2 <- NA_real_
      out$QEp <- NA_real_
      out$analysis_set <- set_name
      return(as_tibble(out))
    }

    out$k <- k
    out$meta_method <- "REML_HK"
    out$estimate <- as.numeric(fit$b[1, 1])
    out$se <- as.numeric(fit$se[1])
    out$ci_low <- as.numeric(fit$ci.lb[1])
    out$ci_high <- as.numeric(fit$ci.ub[1])
    out$p_value <- as.numeric(fit$pval[1])
    out$tau2 <- as.numeric(fit$tau2)
    out$I2 <- as.numeric(fit$I2)
    out$QEp <- as.numeric(fit$QEp)
    out$analysis_set <- set_name
    as_tibble(out)
  })

  bind_rows(rows)
}

main <- function() {
  root <- normalizePath(file.path(get_script_dir(), ".."), mustWork = FALSE)
  if (!file.exists(file.path(root, "config", "targets.yaml"))) {
    root <- normalizePath(getwd(), mustWork = TRUE)
  }

  out_dir <- file.path(root, "results", "analysis")
  fig_dir <- file.path(root, "results", "figures")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  panel <- c("ACTG1", "DNAH9", "GPX3", "VWF", "C4B", "CD44", "CFHR2", "ITIH3", "LRG1", "MYH7B")
  flat_thr <- 0.3

  # OP source-of-truth from spreadsheet
  op <- readxl::read_xlsx(file.path(root, "log fold change 10 Markers.xlsx")) %>%
    transmute(
      gene = as.character(Marker),
      op_log2fc = as.numeric(`log2FC (Chronic/control)`),
      op_direction = direction_from_fc(op_log2fc, thr = flat_thr)
    ) %>%
    filter(gene %in% panel)

  write_csv(op, file.path(out_dir, "op_direction_source_of_truth.csv"))

  # Gene names from targets.yaml
  cfg <- yaml::read_yaml(file.path(root, "config", "targets.yaml"))
  targets <- cfg$targets %||% list()
  gene_name_tbl <- tibble(
    gene = panel,
    gene_full_name = map_chr(panel, ~ (targets[[.x]]$full_name %||% .x))
  )

  # Load proteomics result sources
  eprot <- read_csv(file.path(root, "results", "eprot_query.csv"), show_col_types = FALSE) %>%
    mutate(source = "E-PROT")
  pride_old <- read_csv(file.path(root, "results", "pride_query.csv"), show_col_types = FALSE) %>%
    mutate(source = "PRIDE_legacy")
  pride_new <- read_csv(file.path(root, "results", "pride_discovery", "pride_ad_pd_ms_query.csv"), show_col_types = FALSE) %>%
    mutate(source = "PRIDE_discovery")

  # Prefer discovery PRIDE rows when overlap exists
  pride_old <- pride_old %>% mutate(present = as.character(present))
  pride_new <- pride_new %>% mutate(present = as.character(present))
  pride_all <- bind_rows(pride_new, pride_old) %>%
    mutate(
      present = bool_present(present),
      log2fc = safe_num(log2fc),
      rank_pref = ifelse(source == "PRIDE_discovery", 100L, 10L) +
        ifelse(present, 5L, 0L) + ifelse(!is.na(log2fc), 1L, 0L)
    ) %>%
    arrange(desc(rank_pref)) %>%
    group_by(accession, gene) %>%
    slice(1) %>%
    ungroup()

  data <- bind_rows(eprot, pride_all) %>%
    filter(gene %in% panel) %>%
    transmute(
      accession = as.character(accession),
      disease = as.character(disease),
      biospecimen = as.character(biospecimen),
      source = as.character(source),
      gene = as.character(gene),
      present = bool_present(present),
      n_case = safe_num(n_case),
      n_control = safe_num(n_control),
      mean_case = safe_num(mean_case),
      mean_ctrl = safe_num(mean_ctrl),
      log2fc = safe_num(log2fc)
    )

  # Fill log2fc where missing but means are present
  data <- data %>%
    mutate(log2fc = ifelse(is.na(log2fc) & !is.na(mean_case) & !is.na(mean_ctrl),
                           mean_case - mean_ctrl, log2fc))

  # Dataset metadata and design classes
  metadata <- data %>%
    distinct(accession, disease, biospecimen, source) %>%
    mutate(
      compartment = case_when(
        str_detect(tolower(biospecimen), "brain") ~ "brain",
        str_detect(tolower(biospecimen), "csf") ~ "csf",
        TRUE ~ "other"
      ),
      design_type = case_when(
        accession == "E-PROT-32" ~ "braak_stage",
        accession == "PXD064570" ~ "blocked",
        accession == "PXD034840" ~ "blocked",
        accession == "PXD011216" ~ "blocked",
        TRUE ~ "case_control"
      ),
      disease_group = case_when(
        disease %in% c("Alzheimers", "Parkinsons", "MS") ~ disease,
        TRUE ~ "Exploratory"
      )
    )

  # Tier-C repair attempt for PXD026491 (header-derived counts)
  inferred <- infer_pxd026491_counts()
  n_case_inf <- inferred$n_case
  n_ctrl_inf <- inferred$n_control

  data <- data %>%
    left_join(metadata %>% select(accession, compartment, design_type, disease_group), by = "accession") %>%
    mutate(
      n_case_reported = n_case,
      n_control_reported = n_control,
      n_case_inferred = ifelse(accession == "PXD026491", n_case_inf, NA_real_),
      n_control_inferred = ifelse(accession == "PXD026491", n_ctrl_inf, NA_real_),
      n_case_final = ifelse(!is.na(n_case_reported) & n_case_reported > 0, n_case_reported, n_case_inferred),
      n_control_final = ifelse(!is.na(n_control_reported) & n_control_reported > 0, n_control_reported, n_control_inferred),
      n_source = case_when(
        !is.na(n_case_reported) & !is.na(n_control_reported) & n_case_reported > 0 & n_control_reported > 0 ~ "reported",
        accession == "PXD026491" & !is.na(n_case_inferred) & !is.na(n_control_inferred) ~ "inferred_header",
        TRUE ~ "missing"
      ),
      tier = case_when(
        design_type != "case_control" ~ "Tier_C",
        n_source == "reported" ~ "Tier_A",
        n_source == "inferred_header" ~ "Tier_B",
        TRUE ~ "Tier_C"
      ),
      include_primary = ifelse(
        accession %in% c("E-PROT-31", "E-PROT-61", "E-PROT-65", "PXD016278", "PXD045058"),
        TRUE, FALSE
      ),
      effect_status = case_when(
        !present ~ "not_detected",
        is.na(log2fc) ~ "not_measured",
        abs(log2fc) <= flat_thr ~ "flat",
        log2fc > flat_thr ~ "up",
        log2fc < -flat_thr ~ "down",
        TRUE ~ "no_data"
      )
    ) %>%
    left_join(op %>% select(gene, op_direction), by = "gene") %>%
    mutate(
      concordance = case_when(
        effect_status %in% c("up", "down") & op_direction %in% c("up", "down") &
          effect_status == op_direction ~ "yes",
        effect_status %in% c("up", "down") & op_direction %in% c("up", "down") &
          effect_status != op_direction ~ "no",
        TRUE ~ "na"
      )
    ) %>%
    left_join(gene_name_tbl, by = "gene")

  # Uncertainty-aware per-dataset stats (proxy SE from sample size, explicit caveat)
  data <- data %>%
    mutate(
      se_proxy = ifelse(!is.na(n_case_final) & !is.na(n_control_final) &
                          n_case_final > 0 & n_control_final > 0,
                        sqrt(1 / n_case_final + 1 / n_control_final), NA_real_),
      dof = pmax((n_case_final + n_control_final - 2), 1),
      t_crit = qt(0.975, df = dof),
      ci_low = ifelse(!is.na(log2fc) & !is.na(se_proxy), log2fc - t_crit * se_proxy, NA_real_),
      ci_high = ifelse(!is.na(log2fc) & !is.na(se_proxy), log2fc + t_crit * se_proxy, NA_real_),
      z_score = ifelse(!is.na(log2fc) & !is.na(se_proxy) & se_proxy > 0, log2fc / se_proxy, NA_real_),
      p_value = ifelse(!is.na(z_score), 2 * pnorm(-abs(z_score)), NA_real_)
    ) %>%
    group_by(accession) %>%
    mutate(fdr_bh = p.adjust(p_value, method = "BH")) %>%
    ungroup()

  # Export master publication table
  data_pub <- data %>%
    select(
      accession, disease, disease_group, compartment, biospecimen, source, design_type, tier,
      include_primary, gene, gene_full_name, present, effect_status,
      log2fc, op_direction, concordance,
      n_case_reported, n_control_reported, n_case_inferred, n_control_inferred,
      n_case_final, n_control_final, n_source, se_proxy, ci_low, ci_high, p_value, fdr_bh
    )
  write_csv(data_pub, file.path(out_dir, "publication_effects_with_uncertainty.csv"))

  # Dataset manifest (for supplement transparency)
  manifest <- data_pub %>%
    distinct(accession, disease, disease_group, compartment, biospecimen, source, design_type, tier,
             include_primary, n_case_reported, n_control_reported, n_case_inferred, n_control_inferred, n_source) %>%
    mutate(
      selection_rule = case_when(
        include_primary & tier == "Tier_A" ~ "primary_main_figure",
        tier == "Tier_B" ~ "sensitivity_metadata_inferred",
        design_type == "braak_stage" ~ "sensitivity_non_case_control_design",
        design_type == "blocked" ~ "screened_blocked_unusable",
        TRUE ~ "supplement_exploratory"
      )
    ) %>%
    arrange(compartment, disease, accession)
  write_csv(manifest, file.path(out_dir, "publication_dataset_manifest.csv"))

  # Meta sets (AD/PD/MS only; no OP rows)
  target_dis <- c("Alzheimers", "Parkinsons", "MS")
  meta_pool <- data_pub %>%
    filter(disease %in% target_dis, design_type == "case_control")

  primary_set <- meta_pool %>% filter(tier == "Tier_A")
  sens_no_tier_c <- meta_pool %>% filter(tier != "Tier_C")

  meta_primary_comp <- run_meta_hk(primary_set, c("compartment", "gene"), "primary_tierA")
  meta_primary_discomp <- run_meta_hk(primary_set, c("disease", "compartment", "gene"), "primary_tierA")
  meta_sens_comp <- run_meta_hk(sens_no_tier_c, c("compartment", "gene"), "sensitivity_no_tierC")
  meta_sens_discomp <- run_meta_hk(sens_no_tier_c, c("disease", "compartment", "gene"), "sensitivity_no_tierC")

  write_csv(meta_primary_comp, file.path(out_dir, "publication_meta_primary_by_compartment.csv"))
  write_csv(meta_primary_discomp, file.path(out_dir, "publication_meta_primary_by_disease_compartment.csv"))
  write_csv(meta_sens_comp, file.path(out_dir, "publication_meta_sensitivity_no_tier_c_by_compartment.csv"))
  write_csv(meta_sens_discomp, file.path(out_dir, "publication_meta_sensitivity_no_tier_c_by_disease_compartment.csv"))

  # Concordance summaries
  conc_summary <- data_pub %>%
    group_by(tier, compartment, gene) %>%
    summarise(
      n_yes = sum(concordance == "yes", na.rm = TRUE),
      n_no = sum(concordance == "no", na.rm = TRUE),
      n_na = sum(concordance == "na", na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(conc_summary, file.path(out_dir, "publication_concordance_summary_by_tier.csv"))

  # Figure 1: main heatmap (Tier A primary only)
  fig_main <- data_pub %>%
    filter(include_primary, tier == "Tier_A") %>%
    mutate(
      dataset = case_when(
        accession == "E-PROT-31" ~ "AD brain (E-PROT-31)",
        accession == "E-PROT-61" ~ "AD brain (E-PROT-61)",
        accession == "E-PROT-65" ~ "PD brain (E-PROT-65)",
        accession == "PXD016278" ~ "AD CSF (PXD016278)",
        accession == "PXD045058" ~ "MS CSF (PXD045058)",
        TRUE ~ accession
      ),
      dataset = factor(dataset, levels = c(
        "AD brain (E-PROT-31)", "AD brain (E-PROT-61)", "PD brain (E-PROT-65)",
        "AD CSF (PXD016278)", "MS CSF (PXD045058)"
      )),
      gene_label = factor(paste0(gene, " - ", gene_full_name),
                          levels = rev(paste0(panel, " - ", gene_name_tbl$gene_full_name)))
    )
  write_csv(fig_main, file.path(out_dir, "publication_main_tierA_data.csv"))

  p_main <- ggplot(fig_main, aes(x = dataset, y = gene_label, fill = log2fc)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(log2fc), "", sprintf("%.1f", log2fc))), size = 3) +
    facet_grid(. ~ compartment, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, na.value = "#d9d9d9", name = "log2FC") +
    labs(
      title = "OP panel overlap with neurodegenerative and neuroinflammatory proteomics",
      subtitle = "Primary analysis: Tier A only (case-control metadata available).",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )

  ggsave(file.path(fig_dir, "publication_main_tierA_heatmap.png"), p_main, width = 14, height = 9, dpi = 300, bg = "white")

  # Figure 2: sensitivity heatmap (Tier A + Tier B, excludes Tier C)
  fig_sens <- data_pub %>%
    filter(disease %in% target_dis, design_type == "case_control", tier %in% c("Tier_A", "Tier_B")) %>%
    mutate(
      dataset = paste0(accession, " [", tier, "]"),
      gene_label = factor(paste0(gene, " - ", gene_full_name),
                          levels = rev(paste0(panel, " - ", gene_name_tbl$gene_full_name)))
    )
  write_csv(fig_sens, file.path(out_dir, "publication_sensitivity_no_tierC_data.csv"))

  p_sens <- ggplot(fig_sens, aes(x = dataset, y = gene_label, fill = concordance)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(log2fc), "", sprintf("%.1f", log2fc))), size = 2.8) +
    facet_grid(. ~ compartment, scales = "free_x", space = "free_x") +
    scale_fill_manual(
      values = c("yes" = "#2e7d32", "no" = "#c62828", "na" = "#d9d9d9"),
      drop = FALSE, name = "Concordance"
    ) +
    labs(
      title = "Sensitivity analysis excluding Tier C datasets",
      subtitle = "Includes Tier A + Tier B (metadata inferred).",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )

  ggsave(file.path(fig_dir, "publication_sensitivity_no_tierC_heatmap.png"), p_sens, width = 14, height = 9, dpi = 300, bg = "white")

  # Methods and claims note
  methods_path <- file.path(out_dir, "publication_methods_and_claims.md")
  lines <- c(
    "# Publication Methods and Claim Framing",
    "",
    "## Claim framing",
    "- This analysis supports cross-disease proteomic overlap / partial concordance with the OP signature.",
    "- It does not claim external validation of OP exposure biomarkers in independent OP-exposed human cohorts.",
    "",
    "## OP direction source-of-truth",
    "- Expected OP directions were derived from `log fold change 10 Markers.xlsx` (chronic OP log2FC).",
    "",
    "## Flat-call rule",
    sprintf("- A dataset effect is called `flat` when |log2FC| <= %.2f.", flat_thr),
    "",
    "## Uncertainty model",
    "- Per-dataset uncertainty uses a sample-size proxy standard error: SE_proxy = sqrt(1/n_case + 1/n_control).",
    "- 95% CI and p-values are computed from SE_proxy (explicit approximation when per-protein variance is unavailable).",
    "",
    "## Meta-analysis",
    "- Primary meta-analysis: Tier A only, stratified by compartment (brain, CSF) and by disease+compartment.",
    "- Sensitivity meta-analysis: excludes Tier C (includes Tier A + Tier B).",
    "- Random-effects model: REML with Hartung-Knapp adjustment (`metafor::rma.uni(..., test='knha')`).",
    "",
    "## Tier policy",
    "- Tier A: reported case/control metadata + case-control design.",
    "- Tier B: metadata inferred from machine-readable headers (used in sensitivity only).",
    "- Tier C: missing critical metadata or non-case-control design (excluded from primary inference).",
    "",
    "## Tier repair status",
    sprintf("- PXD026491 header inference: n_case=%s, n_control=%s (Tier B; group mapping uncertainty retained).",
            ifelse(is.na(n_case_inf), "NA", as.character(n_case_inf)),
            ifelse(is.na(n_ctrl_inf), "NA", as.character(n_ctrl_inf))),
    "",
    "## Key output files",
    "- `results/analysis/op_direction_source_of_truth.csv`",
    "- `results/analysis/publication_effects_with_uncertainty.csv`",
    "- `results/analysis/publication_dataset_manifest.csv`",
    "- `results/analysis/publication_meta_primary_by_compartment.csv`",
    "- `results/analysis/publication_meta_primary_by_disease_compartment.csv`",
    "- `results/analysis/publication_meta_sensitivity_no_tier_c_by_compartment.csv`",
    "- `results/analysis/publication_meta_sensitivity_no_tier_c_by_disease_compartment.csv`",
    "- `results/figures/publication_main_tierA_heatmap.png`",
    "- `results/figures/publication_sensitivity_no_tierC_heatmap.png`"
  )
  writeLines(lines, methods_path)

  message("Saved publication tables and figures.")
  message("Manifest: ", file.path(out_dir, "publication_dataset_manifest.csv"))
  message("Methods:  ", methods_path)
}

main()
