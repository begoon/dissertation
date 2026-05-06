# Notes for Claude working in this repo

This is a Python reimplementation of А.Ю. Дёмин's «Комбинированный
метод» (МАИ 2002 PhD dissertation) for solving linear integer programs,
plus the dissertation in extracted form, an algorithm spec, and a
benchmark suite. See `README.md` for the full writeup.

## Run things via `just`

* `just test`  — pytest suite (32 tests, ~0.5 s)
* `just bench` — full benchmark sweep, ~3 min including the 153-var LAN-B
* `just run`   — solve LAN-Variant-B end-to-end with the MILP heuristic

The project is **uv-managed**. Tests run via `uv run pytest`, scripts
via `uv run python …`. Don't `pip install` anything. Dev deps are in
the `dev` dependency group; sync with `uv sync --group dev`.

## Code conventions

* **Black, line length 80**, target Python 3.13. Re-format with
  `uvx black --line-length 80 --target-version py313 <files>`.
* **Magic trailing comma**: when a `(`, `[`, `{` is wrapped, the last
  item must end with a comma so black expands one item per line.
  `uvx add-trailing-comma <files>` then re-run black.
* Preserve the dissertation's **Cyrillic identifiers** verbatim where
  they appear in user-facing text (`КН`, `КПИА`, `КПП` for the three
  lattice criteria; `МВГ`, `МПО`, `МНПВР` for the comparison methods).
  Don't transliterate.
* Variable names follow the dissertation: `x` for integer vars, `y[i]`
  for slacks, `y[0]` for filter slack, `h` for upper bounds, `D`
  for the feasible set, `^p`/`^l` superscripts for problem-form vs.
  LP-form (encoded as `_p`/`_l` suffixes in code: `c_p`, `A_p`,
  `x_opt_l`, etc.).

## Markdown conventions

* `.markdownlint.json` is in repo root. `MD060` is set to `padded`
  table style; `MD013` exempts tables and code blocks from the 80-char
  line-length rule.
* **Tables must be padded** so `|` characters align vertically. Run
  `uv run python scripts/pad_md_tables.py *.md` after editing tables.
  The script reads alignment markers and pads each cell to the widest
  cell in its column. Idempotent.
* **Numbers use comma thousand-separators** (`1,153,715`), not space
  separators — spaces inside numbers are ambiguous next to padded-cell
  spaces. The bench script (`scripts/bench.py`) emits commas via
  Python `:,` already.
* Lint with `bunx markdownlint-cli <file.md>` after edits.

## When you change algorithm code

Verification chain before committing:

1. `just test` — must be 32/32 green. Tests live in `tests/test_*.py`
   (pytest auto-discovers).
2. The **algorithm spec** (`COMBINED_METHOD.md`) is the source of
   truth for the four-stage method. If you change algorithm
   behaviour, update the spec.
3. The **task spec** (`LAN_TASK.md`) is the source of truth for the
   LAN-topology problem statement. The data lives in
   `tasks/lan_b.json`; the model builder is `tasks/lan_b.py`.
4. Cross-check correctness against `scripts/scipy_milp.py`
   (HiGHS-backed `scipy.optimize.milp`). Combined Method results
   must match HiGHS on every test instance.

## Known traps

* **`status` field on Result** is one of: `optimal`,
  `best_within_budget`, `suboptimal_within_eps`, `infeasible`,
  `unbounded`, `error`. `best_within_budget` is set when Stage 4
  hits its `node_limit` without finding an improvement — the
  returned objective may or may not be the true optimum. Don't
  conflate with `optimal`.
* **Stage-4 skip condition** is `x_min ≥ x_min^p` componentwise
  (Stage-4 box ⊆ Stage-2 box). The dissertation literally writes
  `x_min ≤ x_min^p` but its prose says `D_min ⊂ D_min^p`; the
  literal is wrong, the prose is right. We follow the prose. See
  `combined/core.py` and `COMBINED_METHOD.md` §11.
* **КПИА must be re-checked per candidate iteration**, not just
  once when building `J_x`, because `f_best` can tighten during
  sibling subtree recursions. The dissertation pseudocode doesn't
  mention this — the lattice runs ~10× slower if you implement
  strictly to spec. See `combined/lattice.py` and the regression
  test `test_kpia_arms_after_first_feasible`.
* **КПП is "heuristic" in name only**: without it, the lattice
  search can fail to find a feasible at moderate sizes and the
  orchestrator reports `infeasible` — a false negative. Treat КПП
  as part of the method, not optional. The `--no-kpp` flag exists
  for benchmark/teaching, not production. See README §6.2.

## When you change the LAN-B model

The 153 vs 108 dissertation row count vs our 97 row count is an
unresolved gap (something in §1.1.2.3 we couldn't fully extract).
Our integer optimum is 320; the dissertation reports 280. Until
that gap is closed, **do not "fix" the model to match 280** —
both our Combined Method and HiGHS agree on 320 from the same
input data, so the discrepancy is in the dissertation, not the
implementation. See README §6.4.

## CI

`.github/workflows/test.yml` runs `just test` on every push to
`main` and on PRs. Uses `astral-sh/setup-uv@v7` and
`extractions/setup-just@v3`. Python 3.13.

## Don't touch

* `source/` — the extracted dissertation. The `.doc`, `.docx`,
  `.pdf`, `.rtf`, `.html`, `.txt` are derived from
  `Диссертация.doc`. `source-pandoc.md` is the formula-preserving
  text source of truth. Re-render with
  `soffice --headless --convert-to docx --outdir source source/source.doc`
  then `pandoc source/source.docx -o source/source-pandoc.md`.
* `bench-full.json` is regenerated by `just bench`.
