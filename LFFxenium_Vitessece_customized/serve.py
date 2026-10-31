"""Local preview server for site/ that disables browser caching.
If this script fails due to "command not found: uv", instead run 
bash ./serve.sh

Python's plain http.server lets browsers cache files heuristically, so after re-packaging a
sample the browser can keep using stale zarr metadata (.zarray/.zattrs) and render wrongly.

Usage:  uv run python serve.py [--port 8765] [--directory site]
"""

import argparse
import functools
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer


class NoCacheHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-cache")   # always revalidate (cheap 304s)
        super().end_headers()

    def log_message(self, fmt, *args):   # quiet: zarr probing makes many expected 404s
        pass


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--directory", default="site")
    ap.add_argument("--bind", default="127.0.0.1")
    args = ap.parse_args()
    handler = functools.partial(NoCacheHandler, directory=args.directory)
    print(f"serving {args.directory}/ at http://{args.bind}:{args.port}")
    ThreadingHTTPServer((args.bind, args.port), handler).serve_forever()
