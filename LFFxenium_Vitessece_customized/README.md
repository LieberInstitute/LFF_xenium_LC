# This app was produced interactively through Claude Code (Opus 5 and 5.5) using two test-case sample spatialexperiment objects and image files. THIS CODE HAS NOT BEEN LINE-CHECKED BY A HUMAN FOR BUGS OR VULNERABILITIES!
The entire claude back-and-forth that led to the production of all code within this directory is below. All content other than originating SPE objects (qs2 format) and histology/DAPI images were created by Claude. The entire Claude conversation is provided in Claude_DevelopingCustomVites_for_LC_Xenium_Chatlog.pdf .

# Xenium LC viewer

Static Vitessce site showing Xenium cells over H&E or DAPI, coloured by cluster / transferred
label or by gene (logcounts), with LC domain outlines. Everything under `site/` is plain static
files: copy it to any web server.

## Pipeline (per sample, then once for the site)

| Step | Script | Output |
|---|---|---|
| once | `Rscript 00_build_palettes.R` | `viewer_data/label_palettes.json` (paper palettes, renamed to SPE label values) |
| per sample (QC) | `Rscript 01_inspect_inputs.R <spe.qs2> <morphology_focus_0000.ome.tif> viewer_data/<sample>/inspect` | inventory (stdout) + alignment PNGs incl. H&E-vs-DAPI registration check |
| per sample | `Rscript 02_prep_viewer_data.R <spe.qs2> <morphology_focus_0000.ome.tif> viewer_data` | `viewer_data/<sample>/`: `cells.h5ad`, `domains.geojson`, `he_image.tif`, `sample_meta.json`, `qc_domains.png` |
| per sample | `uv run python 03_package_vitessce.py viewer_data/<sample>` | `site/data/<sample>/` (zarr stores, `matrix_info.json`, `vitessce_config.json`), adds the sample (label + title) to `site/samples.json` |
| after all samples | `uv run python 05_global_zscores.py --workers 8` | pooled per-gene mean/SD across samples (`site/zscore_global_stats.json`) and each sample's `layers/zscore_all`; re-run whenever samples are added |
| once | `uv run python 04_vendor_vitessce.py` | `site/vendor/vitessce-4.0.8/` (pinned bundle + one patch, below) |

SPE inputs may be `.qs2` (read with `qs2::qs_read`) or `.RDS`, chosen by file extension.
Sample folders follow `<brnum>/<brnum>.qs2` + `<brnum>/morphology_focus_0000.ome.tif`, e.g.
`Br1039L/Br1039L.qs2`; the trailing `L` (left hemisphere) is part of `brnum`/`sample_id`.

Python dependencies are pinned through `pyproject.toml` / `uv.lock` (`uv sync`).

## Coordinates

Everything is in microns in the Xenium frame (x right, y down, origin at the top-left of the
DAPI image). Cell centroids are already in microns. The DAPI pyramid carries a 0.2125 um/px
scale and the H&E carries 1 / its `imgData` scale factor, in their OME-NGFF metadata.

## Using the viewer

- Every sidebar panel collapses/expands by clicking its title (remembered per browser); the
  others take up the freed space. Panels have no close (x) buttons.
- **Display options** (custom panel, `site/xenviewer-controls.js`):
  - *LC Region Boundaries*: which of the domain definitions is outlined, or None.
  - *Cell size*: marker radius, 1-15 um.
  - *Multi-gene z-score*: whether multi-gene colouring uses z-scores computed within this
    sample or across all samples.
  - *Gene legend* (only while colouring by gene): one gene = logcounts; two or more genes =
    per-cell mean of the genes' z-scored logcounts. Ranges are computed from the same values
    Vitessce colours by.
  - *Colour by <genes>* / *Clear genes*: return to gene colouring after using cell sets, or
    drop the gene selection.
- **Spatial Layers**: H&E, DAPI and cells, each with an inline Alpha slider (hover for
  "Transparency: N% (alpha a)"). DAPI brightness/contrast controls start collapsed (chevron).
- **Cell Sets panel**: pick the grouping (Banksy cluster, transferred SRT / snRNA-seq labels,
  domain membership). Tick/untick sets to show/hide those cells, in either colouring mode.
- **Gene Expression (1 or Multiple)**: plain click adds/removes a gene.
- Sample picker label: `BrNNNN_<APOE>_<Sex>_<age>`; the title over the viewer spells it out.

## Vitessce patches and page customisation

`04_vendor_vitessce.py` vendors Vitessce 4.0.8 (npm integrity verified) and applies four exact
string patches, refusing to run if any target is not found (e.g. after a version bump):
hide cells outside selected sets (show/hide also works with gene colouring); plain click adds
genes in the multi-select list; image channel controls start collapsed; opacity sliders stay
inline in a narrow sidebar.

`site/sidebar.js` + CSS in `site/index.html` re-lay out the right-hand column with CSS grid
(collapsible panels), add the Alpha labels / transparency hints, and put the cell-set level
checkboxes inline. These rely on Vitessce's DOM (role attributes, MUI/rc-tree classes), so
re-check them when upgrading Vitessce. When changing site JS, bump `ASSET_VERSION` in
`site/index.html` so browsers fetch the new modules.

## Notes

- Domain outlines are drawn into the mask (7 um line, max-pooled pyramid) and rendered filled;
  Vitessce's own outline mode has a fixed world-space width that is far too thick zoomed in.
- Caching: browsers must not reuse stale zarr metadata after a sample is re-packaged (a stale
  `.zarray` makes Vitessce read the wrong shape and silently skip data). Serve `site/data/`
  with `Cache-Control: no-cache` (revalidate every time; unchanged files return cheap 304s).
- `site/` is ~485 MB per sample, almost all DAPI (lossless).
- Vitessce's small inset legend averages per-gene display ranges for multi-gene colouring
  (it would read "0 - 238"), so the Display options panel hides it while colouring by gene.
- React is loaded from esm.sh (pinned 18.2.0); vendor it too if the server must avoid CDNs.
- Local preview: `uv run python serve.py` (same as `http.server`, plus `no-cache` headers),
  then open http://localhost:8765
