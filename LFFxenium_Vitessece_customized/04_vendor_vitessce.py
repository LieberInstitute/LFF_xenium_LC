"""04_vendor_vitessce.py

Download the Vitessce JS bundle from npm into site/vendor/ and apply a few small patches, so
the site is served with a pinned, self-hosted copy. Patches (all in one bundle chunk):

1. hide-unselected-cells: in the spatialBeta spot layer, cells that are not in any selected
   cell set are drawn fully transparent (instead of the theme's default grey), so ticking /
   unticking clusters in Cell Sets shows/hides them, also while coloured by gene (the gene
   colormap replaces RGB but keeps the per-cell alpha). No sets selected = all cells shown.
2. click-to-add-gene: with enableMultiSelect, a plain click adds/removes a gene in the Gene
   List (stock behaviour replaces the selection unless Shift is held).
3. channels-collapsed: image-layer channel controls (e.g. DAPI contrast sliders) start
   collapsed in the layer controller.
4. inline-opacity-sliders: the layer controller only moves opacity sliders onto their own
   full-width row when narrower than 250 px (stock: 450 px, i.e. always, in our sidebar).

Patches are exact string replacements and the script fails if a target is not found exactly
once, so a version bump cannot silently drop one.

Usage:  uv run python 04_vendor_vitessce.py [--site site]
"""

import argparse
import base64
import hashlib
import io
import json
import shutil
import tarfile
import urllib.request
from pathlib import Path

VERSION = "4.0.8"
CHUNK = "index-psiMcK01.js"   # bundle chunk holding the spatialBeta layers in 4.0.8

PATCHES = {
    "hide-unselected-cells": ("""}, wRi = (t, e, n) => (i, { index: r }) => {
  const [o, A, s, a] = t && e && t.get(e[r]) || lA(n);
  return [o, A, s, 255 * (a || 1)];
}""", """}, wRi = (t, e, n) => (i, { index: r }) => {
  const xvSetColor = t && e && t.get(e[r]);
  if (!xvSetColor && t && t.size > 0) return [0, 0, 0, 0]; /* xenviewer patch: hide cells outside selected sets */
  const [o, A, s, a] = xvSetColor || lA(n);
  return [o, A, s, 255 * (a || 1)];
}"""),
    # F = Shift held, a = allowMultiple: treat multi-select tables as if Shift were always held
    "click-to-add-gene": ("(J || l) && (!F && (J || !J && a && B.length > 1) ? h([Z])",
                          "(J || l) && (!(F || a) && (J || !J && a && B.length > 1) ? h([Z])"),
    "channels-collapsed": ("[F, Q] = on(!0), { spatialLayerVisible: p, spatialLayerOpacity: S",
                           "[F, Q] = on(!1), { spatialLayerVisible: p, spatialLayerOpacity: S"),
    "inline-opacity-sliders": ("Vbi = 450, jbe", "Vbi = 250, jbe"),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", type=Path, default=Path("site"))
    args = ap.parse_args()

    meta = json.load(urllib.request.urlopen(f"https://registry.npmjs.org/vitessce/{VERSION}"))
    tgz = urllib.request.urlopen(meta["dist"]["tarball"]).read()
    algo, expected = meta["dist"]["integrity"].split("-", 1)
    assert algo == "sha512" and base64.b64encode(hashlib.sha512(tgz).digest()).decode() == expected, \
        "npm integrity check failed"

    out = args.site / "vendor" / f"vitessce-{VERSION}"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    with tarfile.open(fileobj=io.BytesIO(tgz)) as tar:
        for m in tar.getmembers():
            if m.isfile() and m.name.startswith("package/dist/") and m.name.endswith(".js"):
                (out / Path(m.name).name).write_bytes(tar.extractfile(m).read())
            elif m.name in ("package/LICENSE", "package/package.json"):
                (out / Path(m.name).name).write_bytes(tar.extractfile(m).read())

    chunk = out / CHUNK
    src = chunk.read_text()
    for name, (old, new) in PATCHES.items():
        assert src.count(old) == 1, f"patch '{name}' target not found once in {CHUNK}; re-derive it for this version"
        src = src.replace(old, new)
    chunk.write_text(src)
    (out / "XENVIEWER_PATCH.md").write_text(
        f"vitessce {VERSION} from npm (integrity verified), patched in {CHUNK}: "
        + ", ".join(PATCHES) + ".\nApplied by 04_vendor_vitessce.py; see that script for details.\n")
    print(f"vendored vitessce {VERSION} -> {out} (patches: {', '.join(PATCHES)})")


if __name__ == "__main__":
    main()
