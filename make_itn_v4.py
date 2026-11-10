import os, pandas as pd
d = os.path.join("inputs", "interventions")
t = pd.read_csv(os.path.join(d, "ITN_v3.csv"))
t = t[t.year != 2020].copy()
r = t[t.year == 2014].drop_duplicates("node_id", keep="last").copy()
r["year"] = 2017
r["start_day"] = [2661 if n in (1, 2, 3, 4, 8) else 2800 for n in r.node_id]
t = pd.concat([t, r]).sort_values(["node_id", "start_day"])
t.to_csv(os.path.join(d, "ITN_v4.csv"), index=False)
print(t[["node_id", "year", "start_day", "coverage"]].to_string(index=False))
