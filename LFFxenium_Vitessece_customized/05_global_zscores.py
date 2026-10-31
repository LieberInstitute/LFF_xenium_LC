"""05_global_zscores.py

Compute per-gene mean and SD of logcounts pooled over every cell of every packaged sample
(site/samples.json), save them to site/zscore_global_stats.json, and rewrite each sample's
layers/zscore_all with them. Run after packaging all samples, and again whenever samples
are added or re-packaged (03_package_vitessce.py reuses the saved stats for new samples, but
the stats themselves only change when this script is re-run).

Usage:  uv run python 05_global_zscores.py [--site site] [--workers 8]
"""

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import zarr

from importlib import import_module

pkg = import_module("03_package_vitessce")   # reuse zscore() / write_matrix()


def gene_sums(sample_dir):
    info = json.loads((sample_dir / "matrix_info.json").read_text())
    X = zarr.open_group(str(sample_dir / "cells.adata.zarr"), mode="r", zarr_format=2)["X"][:].astype(np.float64)
    return info["genes"], X.shape[0], X.sum(axis=0), (X ** 2).sum(axis=0)


def rewrite_zscore_all(args):
    sample_dir, stats = args
    info_path = sample_dir / "matrix_info.json"
    info = json.loads(info_path.read_text())
    path = sample_dir / "cells.adata.zarr"
    root = zarr.open_group(str(path), mode="a", zarr_format=2)
    gi = [stats["genes"].index(g) for g in info["genes"]]
    z = pkg.zscore(root["X"][:], np.asarray(stats["mean"])[gi], np.asarray(stats["sd"])[gi])
    pkg.write_matrix(root["layers"], "zscore_all", z)
    zarr.consolidate_metadata(str(path), zarr_format=2)
    info["zscore_all_samples"] = stats["samples"]
    info_path.write_text(json.dumps(info))
    return sample_dir.name


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", type=Path, default=Path("site"))
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    samples = json.loads((args.site / "samples.json").read_text())
    dirs = [args.site / "data" / s["id"] for s in samples]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        parts = list(ex.map(gene_sums, dirs))

    genes = parts[0][0]
    for g, *_ in parts[1:]:
        assert g == genes, "gene panels/orders differ between samples"
    n = sum(p[1] for p in parts)
    s1 = sum(p[2] for p in parts)
    s2 = sum(p[3] for p in parts)
    mean = s1 / n
    sd = np.sqrt(np.maximum(s2 / n - mean ** 2, 0))
    stats = {"samples": [s["id"] for s in samples], "n_cells": int(n), "genes": genes,
             "mean": mean.tolist(), "sd": sd.tolist()}
    (args.site / "zscore_global_stats.json").write_text(json.dumps(stats))

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for name in ex.map(rewrite_zscore_all, [(d, stats) for d in dirs]):
            print(f"[{name}] zscore_all updated ({len(stats['samples'])} samples, {n:,} cells)")


if __name__ == "__main__":
    main()
