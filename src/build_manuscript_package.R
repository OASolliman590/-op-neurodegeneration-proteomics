#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  flag <- "--file="
  idx <- grep(flag, args)
  if (length(idx) == 0) return(getwd())
  normalizePath(dirname(sub(flag, "", args[idx][1])), mustWork = TRUE)
}

run_rscript <- function(root, script_name) {
  script_path <- file.path(root, "src", script_name)
  if (!file.exists(script_path)) stop(sprintf("Missing script: %s", script_path))
  status <- system2("Rscript", shQuote(script_path))
  if (!identical(status, 0L)) stop(sprintf("Failed running %s", script_name))
}

safe_move_to_archive <- function(root, rel_path, archive_dir) {
  src <- file.path(root, rel_path)
  if (!file.exists(src)) return(NA_character_)
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  dst <- file.path(archive_dir, basename(src))
  ok <- file.rename(src, dst)
  if (!ok) {
    file.copy(src, dst, overwrite = TRUE)
    file.remove(src)
  }
  dst
}

main <- function() {
  root <- normalizePath(file.path(get_script_dir(), ".."), mustWork = FALSE)
  if (!file.exists(file.path(root, "config", "targets.yaml"))) {
    root <- normalizePath(getwd(), mustWork = TRUE)
  }

  # 1) Rebuild canonical analysis + figures
  run_rscript(root, "publication_overlap_analysis.R")
  run_rscript(root, "mature_grand_finale_ggplot.R")
  run_rscript(root, "null_benchmark_random_panels.R")

  # 2) Build manuscript-facing package directory
  manuscript_root <- file.path(root, "manuscript", "current")
  figs_dir <- file.path(manuscript_root, "figures")
  tables_dir <- file.path(manuscript_root, "tables")
  text_dir <- file.path(manuscript_root, "text")
  qa_dir <- file.path(manuscript_root, "qa")
  dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

  figure_files <- c(
    "results/figures/publication_main_tierA_heatmap.png",
    "results/figures/publication_sensitivity_no_tierC_heatmap.png",
    "results/figures/grand_finale_op_neurodegeneration_ggplot.png",
    "results/figures/grand_finale_op_neurodegeneration_ggplot.pdf",
    "results/figures/null_benchmark_random_panels.png"
  )

  table_files <- c(
    "results/analysis/op_direction_source_of_truth.csv",
    "results/analysis/publication_effects_with_uncertainty.csv",
    "results/analysis/publication_dataset_manifest.csv",
    "results/analysis/publication_main_tierA_data.csv",
    "results/analysis/publication_sensitivity_no_tierC_data.csv",
    "results/analysis/publication_meta_primary_by_compartment.csv",
    "results/analysis/publication_meta_primary_by_disease_compartment.csv",
    "results/analysis/publication_meta_sensitivity_no_tier_c_by_compartment.csv",
    "results/analysis/publication_meta_sensitivity_no_tier_c_by_disease_compartment.csv",
    "results/analysis/null_benchmark_summary.csv"
  )

  text_files <- c(
    "results/analysis/publication_methods_and_claims.md",
    "results/analysis/scientific_validity_appraisal.md"
  )

  copy_manifest <- tibble(
    source = character(),
    destination = character(),
    copied = logical()
  )

  copy_group <- function(paths, target_dir) {
    out <- list()
    for (p in paths) {
      src <- file.path(root, p)
      dst <- file.path(target_dir, basename(p))
      ok <- FALSE
      if (file.exists(src)) {
        ok <- file.copy(src, dst, overwrite = TRUE)
      }
      out[[length(out) + 1]] <- tibble(source = src, destination = dst, copied = isTRUE(ok))
    }
    bind_rows(out)
  }

  copy_manifest <- bind_rows(
    copy_manifest,
    copy_group(figure_files, figs_dir),
    copy_group(table_files, tables_dir),
    copy_group(text_files, text_dir)
  )

  # 3) Hard cleanup of legacy/conflicting figure artifacts
  archive_dir <- file.path(root, "results", "figures", "archive_legacy")
  legacy_figs <- c(
    "results/figures/grand_finale_op_neurodegeneration.png",
    "results/figures/cross_modal_panel_heatmap.png",
    "results/presence_heatmap.png",
    "results/direction_heatmap.png",
    "results/analysis/weighted_concordance_heatmap.png"
  )

  archived <- tibble(original = character(), archived_to = character())
  for (rel in legacy_figs) {
    dst <- safe_move_to_archive(root, rel, archive_dir)
    if (!is.na(dst)) {
      archived <- bind_rows(archived, tibble(original = file.path(root, rel), archived_to = dst))
    }
  }
  existing_archived <- character()
  if (dir.exists(archive_dir)) {
    existing_archived <- list.files(archive_dir, full.names = TRUE)
  }
  if (length(existing_archived) > 0) {
    archived <- bind_rows(
      archived,
      tibble(original = NA_character_, archived_to = existing_archived)
    ) %>%
      distinct(archived_to, .keep_all = TRUE)
  }

  write_csv(copy_manifest, file.path(qa_dir, "package_copy_manifest.csv"))
  write_csv(archived, file.path(qa_dir, "archived_legacy_files.csv"))

  # 4) Run consistency checks after packaging
  run_rscript(root, "manuscript_consistency_check.R")

  # Copy QA report into manuscript package
  qa_reports <- c(
    "results/analysis/manuscript_consistency_checks.csv",
    "results/analysis/manuscript_consistency_report.md",
    "results/analysis/null_benchmark_random_panels.csv",
    "results/analysis/null_benchmark_panel_observed.csv"
  )
  copy_group(qa_reports, qa_dir)

  # 5) Manuscript-facing README
  readme_path <- file.path(manuscript_root, "README.md")
  lines <- c(
    "# Manuscript Package (Current)",
    "",
    sprintf("Build date: %s", as.character(Sys.Date())),
    "",
    "## Claim framing",
    "- This package supports cross-disease proteomic overlap / partial directional concordance with the OP panel.",
    "- It does not claim definitive external validation in independent OP-exposed human cohorts.",
    "",
    "## Main figure",
    "- `figures/publication_main_tierA_heatmap.png`",
    "",
    "## Supplementary figures",
    "- `figures/publication_sensitivity_no_tierC_heatmap.png`",
    "- `figures/grand_finale_op_neurodegeneration_ggplot.pdf`",
    "- `figures/null_benchmark_random_panels.png`",
    "",
    "## Core tables",
    "- `tables/publication_effects_with_uncertainty.csv` (canonical source-of-truth)",
    "- `tables/publication_dataset_manifest.csv`",
    "- `tables/publication_meta_primary_by_compartment.csv`",
    "- `tables/publication_meta_sensitivity_no_tier_c_by_compartment.csv`",
    "- `tables/null_benchmark_summary.csv`",
    "",
    "## QA",
    "- `qa/manuscript_consistency_report.md`",
    "- `qa/manuscript_consistency_checks.csv`",
    "- `qa/archived_legacy_files.csv`"
  )
  writeLines(lines, readme_path)

  message("Manuscript package built at: ", manuscript_root)
  message("Main figure: ", file.path(figs_dir, "publication_main_tierA_heatmap.png"))
}

main()
