# Vitessce customizations

The viewer runs a pinned, self-hosted copy of **Vitessce JS 4.0.8** with a few patches, plus
page-level code that restyles and extends Vitessce's built-in panels. This file lists every
place that depends on Vitessce internals.

> **If Vitessce is updated (JS bundle or the Python `vitessce` package), or any of these files
> is replaced in place, re-check each item below.** Patches target minified code, and the CSS
> and JS target Vitessce's DOM structure. Both can change between releases even when the
> features stay the same. The patch script refuses to run if a target is missing; the
> page-level code fails silently instead (a feature just stops working), so re-test those by
> hand.

## 1. `04_vendor_vitessce.py`: bundle patches

Downloads `vitessce@4.0.8` from npm, checks the npm sha512 integrity hash, extracts `dist/`
into `site/vendor/vitessce-4.0.8/`, and applies four exact-string patches to the chunk
`index-psiMcK01.js`. Each target must occur exactly once or the script stops.

| Patch | What it changes | Why |
|---|---|---|
| `hide-unselected-cells` | Spot colour getter (minified `wRi`, upstream `makeDefaultGetCellColors`): cells not in any selected cell set get alpha 0 instead of the theme's grey. No selection = all cells shown. | Show/hide clusters from Cell Sets, including while colouring by gene (the gene colormap replaces RGB but keeps alpha). |
| `click-to-add-gene` | `SelectableTable` click handler: `!F &&` becomes `!(F \|\| a) &&` (F = Shift held, a = allowMultiple). | With multi-select on, a plain click adds or removes a gene; stock behaviour needs Shift. |
| `channels-collapsed` | Image layer controller expand state `[F, Q] = on(!0)` becomes `on(!1)`. | DAPI brightness/contrast controls start collapsed. |
| `inline-opacity-sliders` | Layer-controller container-query breakpoint `Vbi = 450` becomes `250`. | Below 450 px Vitessce moves opacity sliders to their own full-width row, which is always the case in our sidebar. |

On a version bump, find the equivalent code in the new bundle (the upstream source names are
in the table) and update `CHUNK` and the `PATCHES` strings.

## 2. `site/index.html`

- **Import map** pins React 18.2.0 (esm.sh) and points `vitessce` at the vendored, patched bundle.
- **`ASSET_VERSION`** is appended to the page's own module imports. Bump it whenever
  `sidebar.js` or `xenviewer-controls.js` changes, or browsers keep serving stale copies.
- **Plugin registration**: `pluginViewTypes: [xenviewerControlsPlugin]` on `<Vitessce>`.
- **URL resolution**: data-file URLs and any view prop ending in `Url` are resolved against
  the config's folder, so the site works under any host or sub-path.
- **CSS overrides of Vitessce's DOM**:
  - the sidebar grid layout (`.react-grid-layout` / `.react-grid-item` forced to CSS grid,
    with transforms and transitions disabled)
  - collapse styling, which hooks onto `[role="banner"]`, `[role="heading"]` and `[role="main"]`
  - one-line layer rows (`[class*="layerControllerRoot"] > .MuiGrid-container.MuiGrid-direction-xs-row`,
    MUI `MuiGrid-grid-xs-*` classes)
  - Alpha labels (`.MuiGrid-grid-xs-2:has(.xv-alpha-slider)`)
  - inline cell-set checkboxes (rc-tree `level-0-treenode`, `[class*="levelButtonsContainer"]`,
    `[class*="levelRadioButton"]`)

## 3. `site/sidebar.js`

- Identifies panels by their title text: "Display options", "Spatial Layers", "Cell Sets",
  "Gene Expression…", "Spatial". Renaming a panel in the config means updating `PANELS` here.
- Click-to-collapse on `[role="banner"] [role="heading"]`; state is kept in localStorage
  (`xenviewer.collapsedPanels`).
- Alpha hint: finds `input[type="range"][aria-label^="Adjust opacity"]` inside MUI sliders
  (`.MuiSlider-root` / `.MuiSlider-thumb`) and reads `aria-valuenow`.

## 4. `site/xenviewer-controls.js` (custom "Display options" panel)

- A Vitessce plugin view (`PluginViewType` + `useCoordination`, both exported by `vitessce`).
  It reads and writes these coordination types: `spatialSpotRadius`, `spatialTargetC`,
  `spatialChannelVisible`, `obsType`, `featureSelection`, `obsColorEncoding`,
  `featureValueType`, `featureValueColormapRange`, `legendVisible`.
- Reads expression chunks straight from `cells.adata.zarr` (zlib-compressed, chunk key
  `0.<n>`, layout in `matrix_info.json`) to compute legend ranges. This depends on how
  `03_package_vitessce.py` writes the matrices.
- Hides Vitessce's inset spot legend while colouring by gene, because in 4.0.8 it mislabels
  multi-gene means as "0 – 238".

## 5. `03_package_vitessce.py`: config written for these features

- Uses the Python `vitessce` package (3.9.7, pinned in `uv.lock`) only to build the config JSON.
- View props: `closeButtonVisible=False` on all Vitessce views; Gene list
  `enableMultiSelect=True` and `title="Gene Expression (1 or Multiple)"`; plugin view
  `xenviewerControls` with `domainColumns`, `matrixInfoUrl` and `cellsUrl`.
- Spot layer: `featureAggregationStrategy: "mean"` plus a shared `legendVisible` scope.
  `featureValueType` switches between `expression` (X), `zscore_sample` and `zscore_all`.
- Domain outlines: a single segmentation channel, linked to the spatial view only, so it
  stays out of the layer panel. The control panel switches its `spatialTargetC`.
- Expression matrices are written with zlib so the browser can also decode them (see 4).
- `pyproject.toml` pins `anndata<0.14`, because 0.14 drops Zarr v2 writing, which Vitessce reads.

## Quick re-test checklist after any update

1. Tick and untick Banksy clusters, both while coloured by cluster and by gene: unticked cells
   should disappear.
2. Plain-click two genes: the legend should read "Mean z-score of 2 genes"; flip the z-score toggle.
3. DAPI's channel controls start collapsed and expand cleanly on one line.
4. Opacity sliders are inline with an "Alpha" label, and hovering shows the transparency hint.
5. Each panel collapses and expands from its title, and the others fill the space.
6. LC Region Boundaries dropdown: each option, plus None.
