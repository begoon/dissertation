default: run

run:
    uv run python solve.py tasks/lan_b_zerocost.py --milp-heuristic --node-limit 200000
