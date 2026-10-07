#!/usr/bin/env python
"""Download (optional) and extract text from a PDF.

Usage:
    python tools/pdf2txt.py <pdf-url-or-path> [out.txt] [--max-pages N]

Uses the bundled Python of DSH if run with it:
    C:\\Users\\panyi\\.dsh\\dsh-runtimes\\dsh-primary-runtime\\dependencies\\python\\python.exe
"""
from __future__ import annotations

import sys
import os
import re
import io
import urllib.request

UA = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
)


def fetch(src: str) -> bytes:
    if re.match(r"^https?://", src):
        req = urllib.request.Request(src, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=90) as r:
            return r.read()
    with open(src, "rb") as f:
        return f.read()


def extract(data: bytes, max_pages: int | None = None) -> str:
    from pypdf import PdfReader

    reader = PdfReader(io.BytesIO(data))
    out = []
    for i, page in enumerate(reader.pages):
        if max_pages is not None and i >= max_pages:
            break
        try:
            out.append(f"\n===== page {i + 1} =====\n" + (page.extract_text() or ""))
        except Exception as exc:  # noqa: BLE001
            out.append(f"\n===== page {i + 1} =====\n[extract failed: {exc}]")
    return "".join(out)


_FIX = {
    "\ufb00": "ff", "\ufb01": "fi", "\ufb02": "fl", "\ufb03": "ffi",
    "\ufb04": "ffl", "\ufb05": "st", "\ufb06": "st",
    "\u2018": "'", "\u2019": "'", "\u201c": '"', "\u201d": '"',
    "\u2013": "-", "\u2014": "--", "\u2212": "-", "\u00a0": " ",
    "\u2028": "\n", "\u2029": "\n",
}


def fix_text(text: str) -> str:
    """Make extracted text console/UTF-8 friendly (ligatures, fancy quotes)."""
    for k, v in _FIX.items():
        text = text.replace(k, v)
    return text


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    src = argv[1]
    max_pages = None
    if "--max-pages" in argv:
        max_pages = int(argv[argv.index("--max-pages") + 1])
    out_path = None
    rest = [a for a in argv[2:] if not a.startswith("--") and not a.isdigit()]
    if rest:
        out_path = rest[0]
    text = extract(fetch(src), max_pages)
    text = fix_text(text)
    if out_path:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(text)
        print(f"wrote {out_path} ({len(text)} chars)")
    else:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
