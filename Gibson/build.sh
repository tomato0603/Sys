#!/bin/zsh
# Assemble Gibson.html from the parts.
#
# This is Primer Builder's build.sh, rewritten for Gibson on 2026-09-28. The
# build() function and its three checks are unchanged; only the part list
# and the markers differ.
#
# Why this exists: the engine (parts/shared/engine.js) was once pasted into
# each page by hand, the copies drifted apart, and one release shipped with a
# part silently missing. The concatenation is not trusted on its own: after
# building, the trailing newlines, the line counts and the per-part marker
# strings are checked, and the build fails if any of them is wrong.
#
# Usage:  ./build.sh

set -e
cd "$(dirname "$0")"

fail() { echo "ビルド失敗: $1" >&2; exit 1; }

# Build one page.  $1 = output file, $2... = part paths.
# Each part must carry its marker, listed in MARKERS below under the same path.
build() {
  local out="$1"; shift
  local parts=("$@")
  local tmp="${out}.building"
  local expected=0 part lines

  for part in "${parts[@]}"; do
    [ -f "$part" ] || fail "部品がありません: $part"
    # A part whose last line has no newline glues itself onto the next part's
    # first line. The line-count check below cannot see this (both the part and
    # the output lose the same one line), so it is checked here directly.
    if [ "$(tail -c 1 "$part" | od -An -c | tr -d ' ')" != '\n' ]; then
      fail "$part の末尾に改行がありません（次の部品と 1 行くっつきます）"
    fi
    lines=$(wc -l < "$part")
    expected=$((expected + lines))
  done

  cat "${parts[@]}" > "$tmp"

  # Check 2: no part dropped out of the concatenation, and nothing was written
  # twice. (Check 1 is the trailing-newline test above.)
  local got
  got=$(wc -l < "$tmp")
  if [ "$got" -ne "$expected" ]; then
    rm -f "$tmp"
    fail "$out の行数が合いません（部品の合計 $expected 行、成果物 $got 行）"
  fi

  # Check 3: every part left its marker in the output.
  for part in "${parts[@]}"; do
    local marker="${MARKERS[$part]}"
    [ -n "$marker" ] || fail "目印が定義されていません: $part"
    if ! grep -qF -- "$marker" "$tmp"; then
      rm -f "$tmp"
      fail "$out に $part の目印がありません: $marker"
    fi
  done

  mv "$tmp" "$out"
  echo "$out を書き出しました（$got 行 / 部品 ${#parts[@]} 個）"
}

typeset -A MARKERS
MARKERS=(
  "parts/gibson/head.html"        "<title>Gibson</title>"
  "parts/shared/engine.js"        "/* Primer Builder core -- inverse PCR primer design."
  "gibson.js"                     "var GIBSON = (function () {"
  "parts/shared/test-harness.js"  "var PBTESTS = (function () {"
  "parts/shared/engine-tests.js"  "PBTESTS.EXAMPLE_GENBANK = EXAMPLE_GENBANK;"
  "gibson-tests.js"               "Gibson: 2 本の Tm 差を出す"
  "parts/gibson/ui.js"            "var GIBSONUI = (function () {"
  "parts/gibson/ui-tests.js"      "リクエスト: 何件でも番号順につなげる"
  "parts/gibson/tail.html"        "</html>"
)

if [ -n "$1" ]; then
  fail "引数は取りません（出力は Gibson.html だけです）"
fi

# The order matters: PB before GIBSON (GIBSON calls it), the harness before
# any tests, and ui.js before ui-tests.js because the tests read GIBSONUI.
# engine.js is the same file Primer Builder uses; it is concatenated, never
# copied into another part, so the Tm formula exists in one file only.
build "Gibson.html" \
  parts/gibson/head.html \
  parts/shared/engine.js \
  gibson.js \
  parts/shared/test-harness.js \
  parts/shared/engine-tests.js \
  gibson-tests.js \
  parts/gibson/ui.js \
  parts/gibson/ui-tests.js \
  parts/gibson/tail.html
