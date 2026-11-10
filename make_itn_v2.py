import sys, os, pandas as pd
kill, block, decay, out = float(sys.argv[1]), float(sys.argv[2]), int(float(sys.argv[3])), sys.argv[4]
d = os.path.join("inputs", "interventions")
t = pd.read_csv(os.path.join(d, "ITN.csv"))
t["kill_effect"] = f"PLACEHOLDER_{kill}"; t["block_effect"] = f"PLACEHOLDER_{block}"
t["kill_decay"] = f"PLACEHOLDER_{decay}"; t["block_decay"] = f"PLACEHOLDER_{decay}"
t.to_csv(os.path.join(d, out), index=False)
print("wrote", os.path.join(d, out), "rows", len(t))
print(t[["node_id", "year", "coverage", "kill_effect", "block_effect", "kill_decay", "block_decay"]].head(4).to_string(index=False))
