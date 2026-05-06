default: run

# Run the LAN-B (Variant B) demo end-to-end.
run:
    uv run python solve.py tasks/lan_b_zerocost.py --milp-heuristic --node-limit 200000

# Run the pytest suite.
test:
    uv run pytest tests/

# Run the full benchmark sweep (small problems + 153-var LAN-B).
bench:
    uv run python scripts/bench.py --node-limit 500000 --include-large --json bench-full.json
