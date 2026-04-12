"""
run_queries.py
--------------
Run all query scripts then summarize.

Usage:
  python run_queries.py                  # all sources
  python run_queries.py --only geo       # GEO only
  python run_queries.py --only pride     # PRIDE only
  python run_queries.py --only census    # CellxGene Census only
  python run_queries.py --only expression_atlas # EBI Expression Atlas
  python run_queries.py --only opentargets      # Open Targets evidence
  python run_queries.py --only diseases         # Jensen DISEASES evidence
  python run_queries.py --only hpa       # HPA tissue reference only
  python run_queries.py --only proteomicsdb     # PeptideAtlas detectability
  python run_queries.py --skip census    # skip Census (requires install)
  python run_queries.py --summarize-only # just re-run summarize.py
  python run_queries.py --analyze        # summarize + meta outputs
  python run_queries.py --guardrails     # run guardrail QA gate
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

SRC = Path(__file__).resolve().parent

SOURCES = {
    # Case-control validation sources (produce direction + concordance)
    "geo":               "query_geo.py",
    "pride":             "query_pride.py",
    "eprot":             "query_eprot.py",
    "census":            "query_census.py",
    "expression_atlas":  "query_expression_atlas.py",
    "opentargets":       "query_opentargets.py",
    "diseases":          "query_diseases.py",
    # Reference / detectability layers (produce presence + mean expression only)
    "hpa":               "query_hpa.py",
    "proteomicsdb":      "query_proteomicsdb.py",
}


def run(script: str, extra: list[str] = None):
    path = SRC / script
    cmd  = [sys.executable, str(path)] + (extra or [])
    log.info(f"\n{'─'*50}\n{script}\n{'─'*50}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        log.warning(f"{script} finished with exit code {result.returncode}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--only",  choices=list(SOURCES.keys()),
                        help="Run only this source")
    parser.add_argument("--skip",  nargs="*", choices=list(SOURCES.keys()), default=[],
                        help="Sources to skip")
    parser.add_argument("--summarize-only", action="store_true",
                        help="Skip queries; just re-run summarize.py")
    parser.add_argument("--analyze", action="store_true",
                        help="Run analyze_meta.py after summarize")
    parser.add_argument("--guardrails", action="store_true",
                        help="Run guardrail_check.py after summarize/analyze")
    args = parser.parse_args()

    if args.summarize_only:
        run("summarize.py")
        if args.analyze:
            run("analyze_meta.py")
        if args.guardrails:
            run("guardrail_check.py", ["--strict"])
        return

    sources = ([args.only] if args.only
               else [s for s in SOURCES if s not in args.skip])

    for source in sources:
        run(SOURCES[source])

    # Always summarize at the end
    run("summarize.py")
    if args.analyze:
        run("analyze_meta.py")
    if args.guardrails:
        run("guardrail_check.py", ["--strict"])
    log.info("\nDone. Check results/presence_heatmap.png and results/direction_heatmap.png")


if __name__ == "__main__":
    main()
