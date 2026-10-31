#!/usr/bin/env bash
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

uv run python serve.py --port 8765 --directory site
