#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  flag <- "--file="
  idx <- grep(flag, args)
  if (length(idx) == 0) return(getwd())
  normalizePath(dirname(sub(flag, "", args[idx][1])), mustWork = TRUE)
}

recalc_concordance <- function(effect_status, op_direction) {
  case_when(
    effect_status %in% c("up", "down") & op_direction %in% c("up", "down") & effect_status == op_direction ~ "yes",
    effect_status %in% c("up", "down") & op_direction %in% c("up", "down") & effect_status != op_direction ~ "no",
    TRUE ~ "na"
  )
}

as_status <- function(ok) ifelse(ok, "PASS", "FAIL")
safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

main <- function() {
  root <- normalizePath(file.path(get_script_dir(), ".."), mustWork = FALSE)
  if (!file.exists(file.path(root, "results", "analysis", "publication_effects_with_uncertainty.csv"))) {
    root <- normalizePath(getwd(), mustWork = TRUE)
  }

  analysis_dir <- file.path(root, "results", "analysis")
  manuscript_dir <- file.path(root, "manuscript", "current")
  report_csv <- file.path(analysis_dir, "manuscript_consistency_checks.csv")
  report_md <- file.path(analysis_dir, "manuscript_consistency_report.md")

  canonical_path <- file.path(analysis_dir, "publication_effects_with_uncertainty.csv")
  op_path <- file.path(analysis_dir, "op_direction_source_of_truth.csv")
  req_paths <- c(canonical_path, op_path)
  miss <- req_paths[!file.exists(req_paths)]
  if (length(miss) > 0) stop(sprintf("Missing required file(s): %s", paste(miss, collapse = ", ")))

  canonical <- read_csv(canonical_path, show_col_types = FALSE)
  op <- read_csv(op_path, show_col_types = FALSE) %>% distinct(gene, op_direction)

  checks <- list()

  # 1) Canonical sample-size consistency per accession
  sample_cons <- canonical %>%
    group_by(accession) %>%
    summarise(
      n_case_values = n_distinct(n_case_final[!is.na(n_case_final)]),
      n_ctrl_values = n_distinct(n_control_final[!is.na(n_control_final)]),
      .groups = "drop"
    )
  bad_sample <- sample_cons %>% filter(n_case_values > 1 | n_ctrl_values > 1)
  checks[[length(checks) + 1]] <- tibble(
    check_id = "canonical_sample_size_consistency",
    status = as_status(nrow(bad_sample) == 0),
    details = ifelse(nrow(bad_sample) == 0,
                     "n_case_final/n_control_final are internally consistent per accession.",
                     paste("Mismatched accessions:", paste(bad_sample$accession, collapse = ", ")))
  )

  # 2) Concordance recomputation consistency
  c2 <- canonical %>% mutate(calc = recalc_concordance(effect_status, op_direction))
  bad_c2 <- c2 %>% filter(coalesce(concordance, "na") != calc)
  checks[[length(checks) + 1]] <- tibble(
    check_id = "canonical_concordance_recompute",
    status = as_status(nrow(bad_c2) == 0),
    details = ifelse(nrow(bad_c2) == 0,
                     "Concordance column matches deterministic recomputation.",
                     paste("Rows with mismatch:", nrow(bad_c2)))
  )

  # 3) OP direction source-of-truth lock
  c3 <- canonical %>% distinct(gene, op_direction) %>% left_join(op, by = "gene", suffix = c("_canonical", "_source"))
  bad_c3 <- c3 %>% filter(op_direction_canonical != op_direction_source)
  checks[[length(checks) + 1]] <- tibble(
    check_id = "op_direction_source_lock",
    status = as_status(nrow(bad_c3) == 0),
    details = ifelse(nrow(bad_c3) == 0,
                     "Canonical OP directions match spreadsheet-derived source table.",
                     paste("Direction mismatches in genes:", paste(bad_c3$gene, collapse = ", ")))
  )

  # 4) Cross-artifact count consistency (supplementary figure table)
  ggd_path <- file.path(analysis_dir, "grand_finale_ggplot_data.csv")
  if (file.exists(ggd_path)) {
    ggd <- read_csv(ggd_path, show_col_types = FALSE) %>%
      group_by(accession) %>%
      summarise(n_case_gg = safe_max(n_case), n_ctrl_gg = safe_max(n_control), .groups = "drop")
    can <- canonical %>%
      group_by(accession) %>%
      summarise(n_case_can = safe_max(n_case_final), n_ctrl_can = safe_max(n_control_final), .groups = "drop")
    comp <- inner_join(can, ggd, by = "accession") %>%
      mutate(
        n_case_can = ifelse(is.infinite(n_case_can), NA_real_, n_case_can),
        n_ctrl_can = ifelse(is.infinite(n_ctrl_can), NA_real_, n_ctrl_can),
        n_case_gg = ifelse(is.infinite(n_case_gg), NA_real_, n_case_gg),
        n_ctrl_gg = ifelse(is.infinite(n_ctrl_gg), NA_real_, n_ctrl_gg)
      )
    bad_comp <- comp %>%
      filter((!is.na(n_case_can) & !is.na(n_case_gg) & n_case_can != n_case_gg) |
               (!is.na(n_ctrl_can) & !is.na(n_ctrl_gg) & n_ctrl_can != n_ctrl_gg))

    checks[[length(checks) + 1]] <- tibble(
      check_id = "cross_artifact_count_consistency",
      status = as_status(nrow(bad_comp) == 0),
      details = ifelse(nrow(bad_comp) == 0,
                       "Counts match between canonical table and supplementary ggplot data.",
                       paste("Count mismatch in:", paste(bad_comp$accession, collapse = ", ")))
    )
  }

  # 5) Effect-state taxonomy check
  required_states <- c("not_detected", "flat", "up", "down")
  present_states <- sort(unique(canonical$effect_status))
  miss_states <- setdiff(required_states, present_states)
  checks[[length(checks) + 1]] <- tibble(
    check_id = "effect_state_taxonomy",
    status = as_status(length(miss_states) == 0),
    details = ifelse(length(miss_states) == 0,
                     paste("Effect states present:", paste(present_states, collapse = ", ")),
                     paste("Missing expected states:", paste(miss_states, collapse = ", ")))
  )

  # 5b) Effect-state rule logic consistency
  state_bad <- canonical %>%
    mutate(
      present_bool = coalesce(as.logical(present), FALSE),
      expected_state = case_when(
        !present_bool ~ "not_detected",
        is.na(log2fc) ~ "not_measured",
        abs(log2fc) <= 0.3 ~ "flat",
        log2fc > 0.3 ~ "up",
        log2fc < -0.3 ~ "down",
        TRUE ~ "flat"
      )
    ) %>%
    filter(effect_status != expected_state)
  checks[[length(checks) + 1]] <- tibble(
    check_id = "effect_state_logic_consistency",
    status = as_status(nrow(state_bad) == 0),
    details = ifelse(nrow(state_bad) == 0,
                     "Effect-state calls match deterministic rule mapping.",
                     paste("Rows with state mismatch:", nrow(state_bad)))
  )

  # 6) Language guardrail check in manuscript-facing docs
  text_files <- c(
    file.path(analysis_dir, "publication_methods_and_claims.md"),
    file.path(analysis_dir, "scientific_validity_appraisal.md"),
    file.path(manuscript_dir, "README.md")
  )
  text_files <- text_files[file.exists(text_files)]

  bad_lines <- tibble(file = character(), line_num = integer(), line = character())
  if (length(text_files) > 0) {
    banned <- "\\bvalidated\\b|external validation"
    for (f in text_files) {
      lines <- readLines(f, warn = FALSE)
      hit_idx <- which(str_detect(tolower(lines), banned))
      if (length(hit_idx) > 0) {
        cand <- tibble(file = f, line_num = hit_idx, line = lines[hit_idx]) %>%
          filter(!str_detect(tolower(line), "does not claim external validation|does not claim definitive external validation|not claimed: definitive external validation|unless validated"))
        bad_lines <- bind_rows(bad_lines, cand)
      }
    }
  }

  checks[[length(checks) + 1]] <- tibble(
    check_id = "language_guardrails",
    status = as_status(nrow(bad_lines) == 0),
    details = ifelse(nrow(bad_lines) == 0,
                     "No overclaiming terms in manuscript-facing files.",
                     paste("Flagged lines:", nrow(bad_lines)))
  )

  checks_df <- bind_rows(checks)
  write_csv(checks_df, report_csv)

  md <- c(
    "# Manuscript Consistency Report",
    "",
    sprintf("Date: %s", as.character(Sys.Date())),
    "",
    "## Check Results",
    "",
    "| Check | Status | Details |",
    "|---|---|---|"
  )
  for (i in seq_len(nrow(checks_df))) {
    md <- c(md, sprintf("| %s | %s | %s |", checks_df$check_id[i], checks_df$status[i], checks_df$details[i]))
  }

  if (nrow(bad_lines) > 0) {
    md <- c(md, "", "## Language Guardrail Flags", "", "| File | Line | Text |", "|---|---:|---|")
    for (i in seq_len(nrow(bad_lines))) {
      md <- c(md, sprintf("| %s | %d | %s |", bad_lines$file[i], bad_lines$line_num[i], gsub("\\|", "\\\\|", bad_lines$line[i])))
    }
  }

  writeLines(md, report_md)

  fail_n <- sum(checks_df$status == "FAIL")
  message("Consistency checks complete: ", nrow(checks_df) - fail_n, " PASS / ", fail_n, " FAIL")
  message("Report: ", report_md)

  if (fail_n > 0) quit(status = 1)
}

main()
