#!/usr/bin/env bash
# Usage:
#   bash secondary_genes/Code/run_fisher_int_sec_graph.sh
#   FISHER_N_CORES=12 bash secondary_genes/Code/run_fisher_int_sec_graph.sh
# This script runs fisher_int_sec_graph.R and writes logs to secondary_genes/Code/logs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
R_SCRIPT="${SCRIPT_DIR}/fisher_int_sec_graph.R"
LOG_DIR="${SCRIPT_DIR}/logs"
LOG_FILE="${LOG_DIR}/fisher_int_sec_graph_$(date +%Y-%m-%d_%H%M%S).log"

# Allow override from environment; default to 10 workers.
export FISHER_N_CORES="${FISHER_N_CORES:-10}"

mkdir -p "${LOG_DIR}"
cd "${SCRIPT_DIR}"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Error: Rscript not found in PATH." >&2
  exit 1
fi

echo "Repository root: ${REPO_ROOT}"
echo "Working directory: ${SCRIPT_DIR}"
echo "Parallel workers (requested): ${FISHER_N_CORES}"
echo "Log file: ${LOG_FILE}"
echo "Starting at $(date)"

Rscript "${R_SCRIPT}" 2>&1 | tee "${LOG_FILE}"
status=${PIPESTATUS[0]}

if [[ ${status} -eq 0 ]]; then
  echo "Finished successfully at $(date)"
else
  echo "R script failed with exit code ${status} at $(date)" >&2
fi

exit "${status}"
