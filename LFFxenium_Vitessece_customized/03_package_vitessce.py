"""03_package_vitessce.py

Package one sample's R exports (02_prep_viewer_data.R) into browser-ready files for a static
Vitessce site, and write that sample's Vitessce config.

    site/
      index.html                      viewer page (sample picker; written by hand, not here)
      samples.json                    list of packaged samples (updated by this script)
      data/<sample>/
        cells.adata.zarr              AnnData: X = logcounts (dense float32, gene-chunked, zlib),
                                      layers/zscore_sample, layers/zscore_all (per-gene z-scores
                                      within this sample / across all samples), obsm/spatial (um),
                                      label + domain columns as categoricals
        matrix_info.json              gene order + chunk layout, so the page can compute legend
                                      ranges (incl. mean z-score of several genes) itself
        dapi.ome.zarr                 Xenium DAPI focus image, OME-NGFF 0.4 pyramid (um scale)
        he.ome.zarr                   H&E (RGB) from the SPE imgData, OME-NGFF 0.4 pyramid
        domains.ome.zarr              domain bitmasks: one channel per domain column, 1 = in domain
        vitessce_config.json          config; URLs are relative to this file's folder

All layers share one world coordinate system: microns in the Xenium frame (x right, y down,
origin at the top-left of the DAPI image), carried by each image's OME-NGFF scale transform.

Usage:  uv run python 03_package_vitessce.py viewer_data/Br0942 [--site site]

"All samples" z-scores use site/zscore_global_stats.json when it exists (written by
05_global_zscores.py); otherwise they equal the within-sample z-scores until 05 is run.
"""

import argparse
import json
import os
import re
from pathlib import Path

import anndata as ad
import numpy as np
import tifffile
import zarr
from numcodecs import Blosc, Zlib
from PIL import Image, ImageDraw
from vitessce import CoordinationLevel as CL
from vitessce import VitessceConfig

# --- display defaults -------------------------------------------------------
SPOT_RADIUS_UM = 4.0          # cell centroid marker radius (most DAPI nuclei are ~6-10 um across)
DOMAIN_RASTER_UM = 1.0        # pixel size of the domain outline masks
DOMAIN_LINE_UM = 7.0          # outline thickness
# Outlines are drawn into the mask itself (one channel per domain column) and
# rendered "filled": Vitessce's own bitmask outline mode has a fixed world-space width
# (~37.5 um x stroke x px size), which is too thick zoomed in and vanishes zoomed out.
TILE = 512                    # image chunk edge (px)
GENE_CHUNK = 10               # genes per chunk in X (Vitessce loads one gene column at a time)

LABEL_DISPLAY_NAMES = {
    "BanksyCluster_Label": "Banksy cluster",
    "Xfered_SRT_Label": "SRT label (transferred)",
    "Xfered_snRNA_Label": "snRNA-seq label (transferred)",
}
DEFAULT_COLOR_BY = "BanksyCluster_Label"

COMPRESSOR = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)
# Expression matrices use zlib so the page can also read them with the browser's built-in
# DecompressionStream("deflate") (for the legend), besides Vitessce reading them for colouring.
MATRIX_COMPRESSOR = Zlib(level=6)
DOMAIN_SELECTOR_VIEW = "xenviewerControls"   # plugin view defined in site/index.html


def hex_to_rgb(h):
    h = h.lstrip("#")
    return [int(h[i:i + 2], 16) for i in (0, 2, 4)]


# --- OME-NGFF writing -------------------------------------------------------
def create_level(root, i, shape, dtype):
    return root.create_array(str(i), shape=shape, dtype=dtype,
                             chunks=(1, TILE, TILE), compressors=COMPRESSOR,
                             chunk_key_encoding={"name": "v2", "separator": "/"},
                             # write all-empty chunks too: some static hosts answer a missing
                             # file with 403/HTML instead of 404, which breaks fill-value reads
                             config={"write_empty_chunks": True})


def write_ome_zarr(path, levels, um_per_px, channel_names, channel_colors, windows, name, rgb=False):
    """Write a CYX multiscale image (list of arrays, finest first, each 2x coarser) as
    OME-NGFF 0.4 / Zarr v2 with a micron scale transform per level."""
    root = zarr.open_group(str(path), mode="w", zarr_format=2)
    for i, arr in enumerate(levels):
        create_level(root, i, arr.shape, arr.dtype)[:] = arr
    write_ome_attrs(root, len(levels), um_per_px, channel_names, channel_colors, windows, name, rgb)


def write_ome_attrs(root, n_levels, um_per_px, channel_names, channel_colors, windows, name, rgb=False):
    datasets = [{"path": str(i), "coordinateTransformations": [
        {"type": "scale", "scale": [1.0, um_per_px * 2 ** i, um_per_px * 2 ** i]}]} for i in range(n_levels)]
    root.attrs["multiscales"] = [{
        "version": "0.4", "name": name,
        "axes": [{"name": "c", "type": "channel"},
                 {"name": "y", "type": "space", "unit": "micrometer"},
                 {"name": "x", "type": "space", "unit": "micrometer"}],
        "datasets": datasets,
    }]
    root.attrs["omero"] = {
        "name": name,
        "rdefs": {"model": "color" if rgb else "greyscale"},
        "channels": [{"label": n, "color": c, "active": True,
                      "window": {"start": w[0], "end": w[1], "min": w[0], "max": w[1]}}
                     for n, c, w in zip(channel_names, channel_colors, windows)],
    }


def mean_downsample(a):
    """2x2 mean on a CYX array (edges cropped to even size)."""
    c, h, w = a.shape
    a = a[:, : h - h % 2, : w - w % 2].astype(np.float32)
    return ((a[:, 0::2, 0::2] + a[:, 1::2, 0::2] + a[:, 0::2, 1::2] + a[:, 1::2, 1::2]) / 4)


def pyramid(a, reduce, min_size=TILE):
    levels = [a]
    while max(levels[-1].shape[1:]) > min_size:
        levels.append(reduce(levels[-1]))
    return levels


# --- components -------------------------------------------------------------
def zscore(X, mean, sd):
    """Per-gene z-scores; genes with zero spread get 0."""
    sd = np.where(sd > 0, sd, np.inf)
    return ((X - mean) / sd).astype(np.float32)


def write_matrix(group, name, arr):
    """Dense cells x genes float32 matrix in AnnData's zarr encoding, chunked by gene."""
    z = group.create_array(name, shape=arr.shape, dtype="<f4", chunks=(arr.shape[0], GENE_CHUNK),
                           compressors=MATRIX_COMPRESSOR, overwrite=True,
                           chunk_key_encoding={"name": "v2", "separator": "."})
    z[:] = arr
    z.attrs.update({"encoding-type": "array", "encoding-version": "0.2.0"})


def package_cells(sample_dir, out_dir, meta, global_stats):
    src = ad.read_h5ad(sample_dir / "cells.h5ad")
    keep_obs = ["sample_id", *meta["numeric_obs_columns"], *meta["label_columns"], *meta["domain_columns"]]
    X = np.asarray(src.X.todense(), dtype=np.float32)
    adata = ad.AnnData(
        X=X, obs=src.obs[keep_obs].copy(), var=src.var.copy(),
        obsm={"spatial": src.obsm["spatial"].astype(np.float32)},
    )
    path = out_dir / "cells.adata.zarr"
    ad.settings.zarr_write_format = 2
    adata.write_zarr(path)

    genes = list(adata.var_names)
    z_sample = zscore(X, X.mean(axis=0, dtype=np.float64), X.std(axis=0, dtype=np.float64))
    if global_stats is not None:
        missing = set(genes) - set(global_stats["genes"])
        assert not missing, f"genes missing from global z-score stats: {sorted(missing)[:5]}"
        gi = [global_stats["genes"].index(g) for g in genes]
        z_all = zscore(X, np.asarray(global_stats["mean"])[gi], np.asarray(global_stats["sd"])[gi])
        all_basis = global_stats["samples"]
    else:
        z_all, all_basis = z_sample, [meta["sample_id"]]

    root = zarr.open_group(str(path), mode="a", zarr_format=2)
    write_matrix(root, "X", X)                     # re-write X with zlib + gene chunks
    write_matrix(root["layers"], "zscore_sample", z_sample)
    write_matrix(root["layers"], "zscore_all", z_all)
    zarr.consolidate_metadata(str(path), zarr_format=2)

    (out_dir / "matrix_info.json").write_text(json.dumps({
        "n_obs": int(adata.n_obs), "genes": genes, "gene_chunk": GENE_CHUNK, "dtype": "<f4",
        "compressor": "zlib", "chunk_key": "0.{chunk}",
        "matrices": {"expression": "X", "zscore_sample": "layers/zscore_sample", "zscore_all": "layers/zscore_all"},
        "zscore_all_samples": all_basis,
    }))
    return adata


def package_dapi(meta, out_dir):
    with tifffile.TiffFile(meta["dapi"]["file"]) as tif:
        # XV_THREADS (set by build_site.sh) caps decoding threads when samples run in parallel
        workers = int(os.environ.get("XV_THREADS", "0")) or None
        levels = [lvl.asarray(maxworkers=workers)[None] for lvl in tif.series[0].levels]   # CYX, finest first
    lo, hi = np.percentile(levels[3][levels[3] > 0], [0.5, 99.8])
    write_ome_zarr(out_dir / "dapi.ome.zarr", levels, meta["dapi"]["um_per_px"],
                   ["DAPI"], ["FFFFFF"], [(int(lo), int(hi))], "DAPI")
    return [int(lo), int(hi)]


def package_he(sample_dir, meta, out_dir):
    he = tifffile.imread(sample_dir / meta["he"]["file"])           # HxWx3 uint8
    he = np.moveaxis(he[..., :3], -1, 0)
    levels = [lvl.round().astype(np.uint8) for lvl in pyramid(he, mean_downsample)]
    write_ome_zarr(out_dir / "he.ome.zarr", levels, meta["he"]["um_per_px"],
                   ["R", "G", "B"], ["FF0000", "00FF00", "0000FF"], [(0, 255)] * 3, "H&E", rgb=True)


def max_downsample(a):
    """2x2 max on a CYX array: thin outline lines survive at coarse pyramid levels."""
    c, h, w = a.shape
    a = a[:, : h - h % 2, : w - w % 2]
    return np.maximum.reduce([a[:, 0::2, 0::2], a[:, 1::2, 0::2], a[:, 0::2, 1::2], a[:, 1::2, 1::2]])


def package_domains(sample_dir, meta, out_dir):
    """Rasterise the outlines of the domain polygons (um): one mask channel per domain column,
    1 on a DOMAIN_LINE_UM-wide line along every ring (outer boundary and holes).
    Written channel by channel, so memory stays at one full-resolution mask."""
    gj = json.loads((sample_dir / "domains.geojson").read_text())
    W = int(np.ceil(meta["dapi"]["width_px"] * meta["dapi"]["um_per_px"] / DOMAIN_RASTER_UM))
    H = int(np.ceil(meta["dapi"]["height_px"] * meta["dapi"]["um_per_px"] / DOMAIN_RASTER_UM))
    chans = meta["domain_columns"]
    shapes = [(H, W)]
    while max(shapes[-1]) > TILE:
        h, w = shapes[-1]
        shapes.append(((h - h % 2) // 2, (w - w % 2) // 2))
    root = zarr.open_group(str(out_dir / "domains.ome.zarr"), mode="w", zarr_format=2)
    levels = [create_level(root, i, (len(chans), *sh), np.uint32) for i, sh in enumerate(shapes)]
    rings_by_col = {}
    for f in gj["features"]:
        g = f["geometry"]
        polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
        rings_by_col.setdefault(f["properties"]["domain_set"], []).extend(r for rings in polys for r in rings)
    pixels = {}
    line_px = max(1, round(DOMAIN_LINE_UM / DOMAIN_RASTER_UM))
    for c, dc in enumerate(chans):
        img = Image.new("L", (W, H), 0)
        draw = ImageDraw.Draw(img)
        for ring in rings_by_col.get(dc, []):
            # pixel i spans [i, i+1) um/res, so shift by half a pixel to PIL centres
            xy = [(x / DOMAIN_RASTER_UM - 0.5, y / DOMAIN_RASTER_UM - 0.5) for x, y in ring]
            draw.line(xy + xy[:1], fill=1, width=line_px, joint="curve")
        mask = np.asarray(img, dtype=np.uint32)[None]
        pixels[dc] = int(mask.sum())
        for i, z in enumerate(levels):
            if i:
                mask = max_downsample(mask)   # max-pool: thin lines survive at coarse levels
            z[c] = mask[0]
    write_ome_attrs(root, len(levels), DOMAIN_RASTER_UM, chans, ["FFFFFF"] * len(chans),
                    [(0, 1)] * len(chans), "Domains")
    return pixels


# --- Vitessce config ----------------------------------------------------------
def build_config(sample, adata, meta, palettes, dapi_window):
    label_cols, domain_cols = meta["label_columns"], meta["domain_columns"]
    set_names = {**LABEL_DISPLAY_NAMES, **{dc: f"Domain: {dc}" for dc in domain_cols}}
    set_cols = [*label_cols, *domain_cols]

    vc = VitessceConfig(schema_version="1.0.18", name=sample_label(meta),
                        description=f"{sample}: Xenium cells over H&E / DAPI, with LC domain outlines")
    ds = vc.add_dataset(name=sample)
    ds.add_file(file_type="anndata.zarr", url="cells.adata.zarr",
                options={
                    "obsFeatureMatrix": {"path": "X"},
                    "obsSpots": {"path": "obsm/spatial"},
                    "obsSets": [{"name": set_names[c], "path": f"obs/{c}"} for c in set_cols],
                },
                coordination_values={"obsType": "cell", "featureType": "gene", "featureValueType": "expression"})
    for fvt, path in (("zscore_sample", "layers/zscore_sample"), ("zscore_all", "layers/zscore_all")):
        ds.add_file(file_type="obsFeatureMatrix.anndata.zarr", url="cells.adata.zarr", options={"path": path},
                    coordination_values={"obsType": "cell", "featureType": "gene", "featureValueType": fvt})
    ds.add_file(file_type="image.ome-zarr", url="he.ome.zarr", coordination_values={"fileUid": "he"})
    ds.add_file(file_type="image.ome-zarr", url="dapi.ome.zarr", coordination_values={"fileUid": "dapi"})
    ds.add_file(file_type="obsSegmentations.ome-zarr", url="domains.ome.zarr",
                options={"obsTypesFromChannelNames": True}, coordination_values={"fileUid": "domains"})

    spatial = vc.add_view("spatialBeta", dataset=ds, x=0, y=0, w=8, h=12)
    controls = vc.add_view(DOMAIN_SELECTOR_VIEW, dataset=ds, x=8, y=0, w=4, h=3)
    controls.set_props(domainColumns=domain_cols, matrixInfoUrl="matrix_info.json", cellsUrl="cells.adata.zarr/")
    layers = vc.add_view("layerControllerBeta", dataset=ds, x=8, y=3, w=4, h=3)
    obs_sets = vc.add_view("obsSets", dataset=ds, x=8, y=6, w=4, h=3)
    genes = vc.add_view("featureList", dataset=ds, x=8, y=9, w=4, h=3)
    genes.set_props(enableMultiSelect=True, title="Gene Expression (1 or Multiple)")
    # no close (x) buttons: a closed panel cannot be reopened without reloading; the page
    # makes panels collapsible instead (site/index.html)
    for v in (spatial, layers, obs_sets, genes):
        v.set_props(closeButtonVisible=False)

    # Cell colouring state shared by the spot layer, the cell-set manager and the gene list
    obs_type, feat_type, feat_val_type, sel, colors, enc, feat_sel, cmap_range = vc.add_coordination(
        "obsType", "featureType", "featureValueType", "obsSetSelection", "obsSetColor",
        "obsColorEncoding", "featureSelection", "featureValueColormapRange")
    obs_type.set_value("cell")
    feat_type.set_value("gene")
    feat_val_type.set_value("expression")
    default_levels = list(adata.obs[DEFAULT_COLOR_BY].cat.categories)
    sel.set_value([[set_names[DEFAULT_COLOR_BY], lv] for lv in default_levels])
    colors.set_value([
        {"path": [set_names[c], lv], "color": hex_to_rgb(palettes[c][lv])}
        for c in set_cols for lv in adata.obs[c].cat.categories
    ])
    enc.set_value("cellSetSelection")
    feat_sel.set_value(None)
    cmap_range.set_value([0.0, 1.0])
    for v in (obs_sets, genes):
        v.use_coordination(obs_type, feat_type, feat_val_type, sel, colors, enc, feat_sel, cmap_range)

    # State driven by the custom control panel: cell size, gene legend / z-score mode, and which
    # LC domain outline (one segmentation channel whose target channel the dropdown switches)
    radius, seg_c, seg_visible, seg_obs_type, spot_legend = vc.add_coordination(
        "spatialSpotRadius", "spatialTargetC", "spatialChannelVisible", "obsType", "legendVisible")
    radius.set_value(SPOT_RADIUS_UM)
    # Vitessce's inset legend mislabels multi-gene aggregates (it averages per-gene display
    # extents), so the control panel hides it while colouring by gene and shows its own.
    spot_legend.set_value(True)
    seg_c.set_value(0)
    seg_visible.set_value(True)
    seg_obs_type.set_value(domain_cols[0])
    controls.use_coordination(radius, seg_c, seg_visible, seg_obs_type, spot_legend,
                              feat_sel, enc, feat_val_type, cmap_range)

    outline_rgb = hex_to_rgb(palettes["domain_outline_colour"])
    vc.link_views_by_dict([spatial, layers], {
        "imageLayer": CL([
            {
                "fileUid": "he",
                "spatialLayerVisible": True,
                "spatialLayerOpacity": 1.0,
                "photometricInterpretation": "RGB",
                "imageChannel": CL([
                    {"spatialTargetC": i, "spatialChannelColor": rgb, "spatialChannelVisible": True,
                     "spatialChannelOpacity": 1.0, "spatialChannelWindow": [0, 255]}
                    for i, rgb in enumerate([[255, 0, 0], [0, 255, 0], [0, 0, 255]])
                ]),
            },
            {
                "fileUid": "dapi",
                "spatialLayerVisible": False,
                "spatialLayerOpacity": 1.0,
                "photometricInterpretation": "BlackIsZero",
                "imageChannel": CL([
                    {"spatialTargetC": 0, "spatialChannelColor": [255, 255, 255], "spatialChannelVisible": True,
                     "spatialChannelOpacity": 1.0, "spatialChannelWindow": dapi_window},
                ]),
            },
        ]),
        "spotLayer": CL([
            {
                "obsType": obs_type,
                "spatialLayerVisible": True,
                "spatialLayerOpacity": 1.0,
                "spatialSpotRadius": radius,
                "spatialSpotFilled": True,
                "spatialSpotStrokeWidth": 1.0,
                "spatialLayerColor": [255, 255, 255],
                "featureType": feat_type,
                "featureValueType": feat_val_type,
                "obsColorEncoding": enc,
                "featureSelection": feat_sel,
                "featureValueColormap": "viridis",
                "featureValueColormapRange": cmap_range,
                "featureAggregationStrategy": "mean",
                "obsSetSelection": sel,
                "obsSetColor": colors,
                "obsHighlight": None,
                "tooltipsVisible": True,
                "legendVisible": spot_legend,
            },
        ]),
    }, meta=True)

    # Domain outlines: linked to the spatial view only (kept out of the layer panel); the
    # control panel's dropdown sets the channel shown, or hides it.
    vc.link_views_by_dict([spatial], {
        "segmentationLayer": CL([
            {
                "fileUid": "domains",
                "spatialLayerVisible": True,
                "spatialLayerOpacity": 1.0,
                "segmentationChannel": CL([
                    {
                        "obsType": seg_obs_type,
                        "spatialTargetC": seg_c,
                        "spatialChannelVisible": seg_visible,
                        "spatialChannelOpacity": 1.0,
                        "spatialChannelColor": outline_rgb,
                        "obsColorEncoding": "spatialChannelColor",
                        "spatialSegmentationFilled": True,
                        "spatialSegmentationStrokeWidth": 1.0,
                        "obsHighlight": None,
                    },
                ]),
            },
        ]),
    }, meta=True)

    # start centred on the tissue, whole section in view
    x0, x1 = meta["centroid_range_um"]["x"]
    y0, y1 = meta["centroid_range_um"]["y"]
    for ct, val in (("spatialTargetX", (x0 + x1) / 2), ("spatialTargetY", (y0 + y1) / 2), ("spatialZoom", -3.6)):
        scope, = vc.add_coordination(ct)
        scope.set_value(val)
        for v in (spatial, layers):
            v.use_coordination(scope)
    return vc.to_dict(base_url="")


def apoe_long(code):
    """'E4hom' -> 'APOE4 homozygous', 'E3E4' / 'E3/E4' -> 'APOE3/4 heterozygous'."""
    alleles = re.findall(r"[234]", str(code))
    if "hom" in str(code).lower() and len(alleles) == 1:
        return f"APOE{alleles[0]} homozygous"
    if len(alleles) == 2:
        a, b = sorted(alleles)
        return f"APOE{a} homozygous" if a == b else f"APOE{a}/{b} heterozygous"
    return f"APOE {code}"


def sample_label(meta):
    """Dropdown label: BrNNNN_<APOE>_<Sex>_<age>."""
    d = meta["donor"]
    return f"{d['brnum']}_{d['apoe']}_{d['sex']}_{int(np.floor(float(d['age']) + 0.5))}"


def sample_title(meta):
    """Title over the viewer: 'Br0942 (APOE4 homozygous, Male, 62 y/o)'."""
    d = meta["donor"]
    sex = {"M": "Male", "F": "Female"}.get(str(d["sex"]).upper(), str(d["sex"]))
    return f"{d['brnum']} ({apoe_long(d['apoe'])}, {sex}, {int(np.floor(float(d['age']) + 0.5))} y/o)"


def write_index_entry(out_dir, sample, meta):
    """Each sample writes its own entry; samples.json is rebuilt from all entries, so parallel
    packaging runs never overwrite each other's additions."""
    entry = {"id": sample, "label": sample_label(meta), "title": sample_title(meta),
             "config": f"data/{sample}/vitessce_config.json", "n_cells": meta["n_cells"]}
    (out_dir / "sample_index_entry.json").write_text(json.dumps(entry, indent=2))


def rebuild_sample_index(site):
    """site/samples.json = every site/data/<sample>/sample_index_entry.json (plus entries already
    listed for packaged samples that predate entry files), sorted by id; written atomically."""
    idx_path = site / "samples.json"
    entries = {}
    if idx_path.exists():
        for s in json.loads(idx_path.read_text()):
            if (site / "data" / s["id"] / "vitessce_config.json").exists():
                entries[s["id"]] = s
    for p in (site / "data").glob("*/sample_index_entry.json"):
        e = json.loads(p.read_text())
        entries[e["id"]] = e
    tmp = idx_path.with_name(f".samples.json.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(sorted(entries.values(), key=lambda s: s["id"]), indent=2))
    os.replace(tmp, idx_path)
    return len(entries)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample_dir", type=Path, nargs="?", help="viewer_data/<sample> folder written by 02_prep_viewer_data.R")
    ap.add_argument("--site", type=Path, default=Path("site"))
    ap.add_argument("--skip-images", action="store_true", help="reuse existing image stores (config/cells only)")
    ap.add_argument("--rebuild-index", action="store_true",
                    help="only rebuild site/samples.json from the per-sample index entries")
    args = ap.parse_args()

    if args.rebuild_index:
        print(f"site/samples.json: {rebuild_sample_index(args.site)} samples")
        return
    if args.sample_dir is None:
        ap.error("sample_dir is required unless --rebuild-index is given")
    meta = json.loads((args.sample_dir / "sample_meta.json").read_text())
    palettes = json.loads((args.sample_dir / meta["palettes"]).read_text())
    sample = meta["sample_id"]
    out_dir = args.site / "data" / sample
    out_dir.mkdir(parents=True, exist_ok=True)

    stats_path = args.site / "zscore_global_stats.json"
    global_stats = json.loads(stats_path.read_text()) if stats_path.exists() else None
    print(f"[{sample}] cells"); adata = package_cells(args.sample_dir, out_dir, meta, global_stats)
    if args.skip_images:
        dapi_window = json.loads((out_dir / "package_info.json").read_text())["dapi_window"]
    else:
        print(f"[{sample}] H&E"); package_he(args.sample_dir, meta, out_dir)
        print(f"[{sample}] domains"); dom_px = package_domains(args.sample_dir, meta, out_dir)
        print(f"[{sample}] DAPI (decoding JPEG 2000 pyramid)"); dapi_window = package_dapi(meta, out_dir)
        (out_dir / "package_info.json").write_text(json.dumps(
            {"dapi_window": dapi_window, "domain_pixels": dom_px, "domain_raster_um": DOMAIN_RASTER_UM,
             "domain_line_um": DOMAIN_LINE_UM}, indent=2))

    cfg = build_config(sample, adata, meta, palettes, dapi_window)
    (out_dir / "vitessce_config.json").write_text(json.dumps(cfg, indent=2))
    write_index_entry(out_dir, sample, meta)
    rebuild_sample_index(args.site)
    print(f"[{sample}] done -> {out_dir}")


if __name__ == "__main__":
    main()
