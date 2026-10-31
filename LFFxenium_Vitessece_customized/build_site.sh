#!/usr/bin/env bash
# build_site.sh: batch driver for the Xenium LC viewer.
#
# Finds every directory starting with "Br" next to this script, and for each one runs
#   01_inspect_inputs.R    (QC: inventory + H&E/DAPI registration plots)
#   02_prep_viewer_data.R  (R prep: cells, domain outlines, H&E, metadata)
#   03_package_vitessce.py (browser-ready zarr stores + Vitessce config)
# several samples at a time, then the once-per-site steps:
#   rebuild site/samples.json, and 05_global_zscores.py (z-scores pooled across all samples).
# One-time setup (Python env, colour palettes, vendored Vitessce) runs first, only if needed.
#
# Each sample folder needs:
#   <Br...>/morphology_focus_0000.ome.tif
#   <Br...>/<Br...>.qs2   (or <Br...>.RDS, or else the folder's only .qs2/.RDS file)
#
# Usage:  ./build_site.sh [options] [SAMPLE ...]
#   SAMPLE ...          only these sample folders (default: every Br* folder)
#   -j, --jobs N        samples processed at once (default 4)
#   -t, --threads N     threads per sample (default 2; keep jobs x threads <= cores you can spare)
#   --skip-done         skip samples already packaged (site/data/<sample>/vitessce_config.json exists)
#   --skip-qc           skip 01_inspect_inputs.R
#   --rebuild-palettes  re-seed viewer_data/label_palettes.json from the paper palettes
#                       (otherwise kept, including colours auto-assigned for new labels)
#   --revendor          re-download and re-patch the Vitessce bundle (needs internet)
#   -n, --dry-run       list samples and what would run, then stop
#   -h, --help          this help
#
# Logs: logs/<sample>.log per sample (QC inventory: viewer_data/<sample>/inspect/).

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

VITESSCE_VERSION="4.0.8"   # keep in step with VERSION in 04_vendor_vitessce.py
DAPI_NAME="morphology_focus_0000.ome.tif"
LOG_DIR="logs"

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
die()  { echo "ERROR: $*" >&2; exit 1; }
note() { echo "[$(date +%H:%M:%S)] $*"; }

# SPE file for a sample folder: <name>.qs2, <name>.RDS/.rds, else the folder's only .qs2/.RDS
find_spe() {
    local dir="$1" name f
    name="$(basename "$dir")"
    for f in "$dir/$name.qs2" "$dir/$name.RDS" "$dir/$name.rds"; do
        [[ -f "$f" ]] && { echo "$f"; return 0; }
    done
    local found=() n=0
    for f in "$dir"/*.qs2 "$dir"/*.RDS "$dir"/*.rds; do
        [[ -f "$f" ]] && { found+=("$f"); n=$((n + 1)); }
    done
    [[ $n -eq 1 ]] && { echo "${found[0]}"; return 0; }
    return 1
}

# Limit BLAS / OpenMP / blosc / data.table threads so parallel samples don't oversubscribe
limit_threads() {
    export XV_THREADS="$1" OMP_NUM_THREADS="$1" OPENBLAS_NUM_THREADS="$1" MKL_NUM_THREADS="$1" \
           VECLIB_MAXIMUM_THREADS="$1" NUMEXPR_NUM_THREADS="$1" BLOSC_NTHREADS="$1" \
           R_DATATABLE_NUM_THREADS="$1"
}

# ---------------------------------------------------------------------------
# per-sample worker (invoked by xargs as: build_site.sh --one <sample>)
# ---------------------------------------------------------------------------
run_sample() {
    local name="$1" dir="$1" spe dapi log t0 step
    log="$LOG_DIR/$name.log"
    t0=$SECONDS
    limit_threads "${XV_THREADS:-2}"
    spe="$(find_spe "$dir")" || { note "$name: FAILED (no .qs2/.RDS found)"; echo "no SPE file" > "$log"; return 1; }
    dapi="$dir/$DAPI_NAME"
    : > "$log"
    note "$name: started"

    if [[ "${XV_SKIP_QC:-0}" != 1 ]]; then
        step="01 QC"
        mkdir -p "viewer_data/$name/inspect"
        { echo "== $step"; Rscript 01_inspect_inputs.R "$spe" "$dapi" "viewer_data/$name/inspect" \
            > "viewer_data/$name/inspect/01_inspect_inputs.log" 2>&1; } >> "$log" 2>&1 \
            || { note "$name: FAILED at $step (see viewer_data/$name/inspect/01_inspect_inputs.log)"; return 1; }
    fi
    step="02 R prep"
    { echo "== $step"; Rscript 02_prep_viewer_data.R "$spe" "$dapi" viewer_data; } >> "$log" 2>&1 \
        || { note "$name: FAILED at $step (see $log)"; return 1; }
    step="03 packaging"
    { echo "== $step"; uv run --frozen python 03_package_vitessce.py "viewer_data/$name"; } >> "$log" 2>&1 \
        || { note "$name: FAILED at $step (see $log)"; return 1; }

    note "$name: done in $((SECONDS - t0)) s"
}

if [[ "${1:-}" == "--one" ]]; then
    run_sample "$2"
    exit $?
fi

# ---------------------------------------------------------------------------
# options
# ---------------------------------------------------------------------------
JOBS=4; THREADS=2; SKIP_DONE=0; SKIP_QC=0; REBUILD_PAL=0; REVENDOR=0; DRY_RUN=0
REQUESTED=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)          JOBS="${2:?}"; shift 2 ;;
        -t|--threads)       THREADS="${2:?}"; shift 2 ;;
        --skip-done)        SKIP_DONE=1; shift ;;
        --skip-qc)          SKIP_QC=1; shift ;;
        --rebuild-palettes) REBUILD_PAL=1; shift ;;
        --revendor)         REVENDOR=1; shift ;;
        -n|--dry-run)       DRY_RUN=1; shift ;;
        -h|--help)          sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)                 die "unknown option $1 (see --help)" ;;
        *)                  REQUESTED+=("${1%/}"); shift ;;
    esac
done
[[ "$JOBS" =~ ^[1-9][0-9]*$ && "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "--jobs and --threads must be positive integers"

# ---------------------------------------------------------------------------
# discover samples
# ---------------------------------------------------------------------------
CANDIDATES=()
if [[ ${#REQUESTED[@]} -gt 0 ]]; then
    CANDIDATES=("${REQUESTED[@]}")
else
    for d in Br*/; do [[ -d "$d" ]] && CANDIDATES+=("${d%/}"); done
fi
[[ ${#CANDIDATES[@]} -gt 0 ]] || die "no Br* sample folders found in $SCRIPT_DIR"

SAMPLES=()
note "found ${#CANDIDATES[@]} sample folder(s):"
for s in "${CANDIDATES[@]}"; do
    if [[ ! -d "$s" ]]; then echo "  $s: SKIP (no such folder)"; continue; fi
    if ! spe="$(find_spe "$s")"; then echo "  $s: SKIP (no single .qs2/.RDS file)"; continue; fi
    if [[ ! -f "$s/$DAPI_NAME" ]]; then echo "  $s: SKIP (no $DAPI_NAME)"; continue; fi
    if [[ $SKIP_DONE -eq 1 && -f "site/data/$s/vitessce_config.json" ]]; then
        echo "  $s: skip (already packaged; --skip-done)"; continue
    fi
    echo "  $s: $spe + $s/$DAPI_NAME"
    SAMPLES+=("$s")
done

if [[ $DRY_RUN -eq 1 ]]; then
    echo
    echo "dry run: would run, with $JOBS sample(s) at a time x $THREADS thread(s) each:"
    [[ -f viewer_data/label_palettes.json && $REBUILD_PAL -eq 0 ]] || echo "  Rscript 00_build_palettes.R"
    [[ -f "site/vendor/vitessce-$VITESSCE_VERSION/index.min.js" && $REVENDOR -eq 0 ]] || echo "  uv run python 04_vendor_vitessce.py"
    for s in ${SAMPLES[@]+"${SAMPLES[@]}"}; do
        [[ $SKIP_QC -eq 1 ]] || echo "  [$s] Rscript 01_inspect_inputs.R ..."
        echo "  [$s] Rscript 02_prep_viewer_data.R ... ; uv run python 03_package_vitessce.py viewer_data/$s"
    done
    echo "  uv run python 03_package_vitessce.py --rebuild-index"
    echo "  uv run python 05_global_zscores.py"
    exit 0
fi

# ---------------------------------------------------------------------------
# one-time setup
# ---------------------------------------------------------------------------
# Tools: use PATH, else the usual install folders (uv's installer puts it in ~/.local/bin,
# which shells such as conda's don't always have on PATH).
find_tool() {
    command -v "$1" >/dev/null && return 0
    local d
    for d in "$HOME/.local/bin" "$HOME/.cargo/bin" /opt/homebrew/bin /usr/local/bin \
             /Library/Frameworks/R.framework/Resources/bin; do
        [[ -x "$d/$1" ]] && { export PATH="$d:$PATH"; return 0; }
    done
    return 1
}
find_tool Rscript || die "Rscript not found (install R, or add its bin folder to PATH)"
find_tool uv || die "uv not found on PATH or in ~/.local/bin (https://docs.astral.sh/uv/)"
for f in index.html xenviewer-controls.js sidebar.js; do
    [[ -f "site/$f" ]] || die "site/$f missing: copy the viewer page files (site/index.html, site/xenviewer-controls.js, site/sidebar.js) into site/ here; they are hand-written, not generated"
done
mkdir -p "$LOG_DIR" viewer_data site/data
T_START=$SECONDS

note "Python environment (uv sync)"
uv sync --frozen --quiet || die "uv sync failed"

if [[ ! -f viewer_data/label_palettes.json || $REBUILD_PAL -eq 1 ]]; then
    note "colour palettes (00_build_palettes.R)"
    Rscript 00_build_palettes.R > "$LOG_DIR/00_build_palettes.log" 2>&1 \
        || die "00_build_palettes.R failed (see $LOG_DIR/00_build_palettes.log)"
fi
rm -rf viewer_data/label_palettes.json.lock   # stale lock from an interrupted run

if [[ ! -f "site/vendor/vitessce-$VITESSCE_VERSION/index.min.js" || $REVENDOR -eq 1 ]]; then
    note "Vitessce bundle (04_vendor_vitessce.py; downloads from npm)"
    uv run --frozen python 04_vendor_vitessce.py > "$LOG_DIR/04_vendor_vitessce.log" 2>&1 \
        || die "04_vendor_vitessce.py failed (see $LOG_DIR/04_vendor_vitessce.log)"
fi

# ---------------------------------------------------------------------------
# per-sample steps, in parallel
# ---------------------------------------------------------------------------
FAILED=0
if [[ ${#SAMPLES[@]} -gt 0 ]]; then
    note "processing ${#SAMPLES[@]} sample(s): $JOBS at a time, $THREADS thread(s) each"
    export XV_THREADS="$THREADS" XV_SKIP_QC="$SKIP_QC"
    printf '%s\n' "${SAMPLES[@]}" | xargs -n 1 -P "$JOBS" bash "$SCRIPT_DIR/build_site.sh" --one || FAILED=1
else
    note "no samples to process"
fi

# ---------------------------------------------------------------------------
# once-per-site steps (run over everything packaged so far, even if some samples failed)
# ---------------------------------------------------------------------------
note "sample index (site/samples.json)"
uv run --frozen python 03_package_vitessce.py --rebuild-index || die "rebuilding samples.json failed"

if compgen -G "site/data/*/matrix_info.json" > /dev/null; then
    WORKERS=$((JOBS * THREADS)); [[ $WORKERS -gt 8 ]] && WORKERS=8
    note "z-scores across all samples (05_global_zscores.py)"
    uv run --frozen python 05_global_zscores.py --workers "$WORKERS" > "$LOG_DIR/05_global_zscores.log" 2>&1 \
        || die "05_global_zscores.py failed (see $LOG_DIR/05_global_zscores.log)"
    tail -n 1 "$LOG_DIR/05_global_zscores.log"
fi

note "finished in $((SECONDS - T_START)) s"
if [[ $FAILED -ne 0 ]]; then
    echo "Some samples FAILED; see the lines above and $LOG_DIR/<sample>.log" >&2
    exit 1
fi
echo "Preview: uv run python serve.py  ->  http://localhost:8765"
