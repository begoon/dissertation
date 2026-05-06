"""Pad every Markdown table in a file so column pipes align vertically.

Each cell is padded to the width of its widest cell in the column, on the
same side as its alignment marker (`---` left, `---:` right, `:---:` center).
The separator row uses `-` chars padded to match.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

PIPE_LINE = re.compile(r"^\s*\|.*\|\s*$")
SEP_CELL = re.compile(r"^\s*:?-+:?\s*$")


def split_row(line: str) -> list[str]:
    # Drop the leading and trailing pipe, split on `|` (no escapes used here).
    inner = line.strip()
    assert inner.startswith("|") and inner.endswith("|")
    return [c.strip() for c in inner[1:-1].split("|")]


def is_separator(cells: list[str]) -> bool:
    return all(SEP_CELL.match(c or "") for c in cells)


def alignment(sep_cell: str) -> str:
    s = sep_cell.strip()
    left = s.startswith(":")
    right = s.endswith(":")
    if left and right:
        return "center"
    if right:
        return "right"
    return "left"


def visual_len(s: str) -> int:
    # Treat each char as 1 column. Markdown editors render combining marks
    # narrowly anyway; this is good enough for our tables.
    return len(s)


def pad_cell(content: str, width: int, align: str) -> str:
    pad = width - visual_len(content)
    if pad <= 0:
        return content
    if align == "right":
        return " " * pad + content
    if align == "center":
        left = pad // 2
        right = pad - left
        return " " * left + content + " " * right
    return content + " " * pad


def format_table(rows: list[list[str]], aligns: list[str]) -> list[str]:
    n = len(aligns)
    widths = [0] * n
    for row in rows:
        for i, cell in enumerate(row):
            if i < n:
                widths[i] = max(widths[i], visual_len(cell))
    # Separator has at least 3 dashes; widen to match column.
    out: list[str] = []
    for idx, row in enumerate(rows):
        if idx == 1:
            # Separator row: compose `-`s with optional `:` markers.
            parts = []
            for i, a in enumerate(aligns):
                w = max(widths[i], 3)
                if a == "right":
                    parts.append("-" * (w - 1) + ":")
                elif a == "center":
                    parts.append(":" + "-" * (w - 2) + ":")
                else:
                    parts.append("-" * w)
            out.append("| " + " | ".join(parts) + " |")
        else:
            cells = [
                pad_cell(row[i] if i < len(row) else "", widths[i], aligns[i])
                for i in range(n)
            ]
            out.append("| " + " | ".join(cells) + " |")
    return out


def reformat(text: str) -> str:
    lines = text.split("\n")
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if not PIPE_LINE.match(line):
            out.append(line)
            i += 1
            continue
        # Collect all consecutive table rows
        j = i
        rows: list[list[str]] = []
        while j < len(lines) and PIPE_LINE.match(lines[j]):
            rows.append(split_row(lines[j]))
            j += 1
        # Need at least 2 rows (header + separator) and second must be a separator.
        if len(rows) < 2 or not is_separator(rows[1]):
            out.extend(lines[i:j])
            i = j
            continue
        aligns = [alignment(c) for c in rows[1]]
        # Pad shorter rows to the alignment-row column count
        n_cols = len(aligns)
        rows = [r + [""] * (n_cols - len(r)) if len(r) < n_cols else r[:n_cols]
                for r in rows]
        # Indentation prefix taken from the first line (for tables under list items)
        indent_match = re.match(r"^(\s*)", lines[i])
        prefix = indent_match.group(1) if indent_match else ""
        for row_line in format_table(rows, aligns):
            out.append(prefix + row_line)
        i = j
    return "\n".join(out)


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: pad_md_tables.py FILE [FILE ...]", file=sys.stderr)
        return 2
    for path in sys.argv[1:]:
        p = Path(path)
        new = reformat(p.read_text())
        if new != p.read_text():
            p.write_text(new)
            print(f"padded: {path}")
        else:
            print(f"unchanged: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
