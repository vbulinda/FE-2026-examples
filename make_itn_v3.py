import os, pandas as pd
d = os.path.join("inputs", "interventions")
t = pd.read_csv(os.path.join(d, "ITN_v2.csv"))
P1 = [1, 2, 3, 4, 8]
t = t[~(t.node_id.isin(P1) & (t.year == 2015))].copy()
t.loc[t.node_id.isin(P1) & (t.year == 2014), "start_day"] = 1730
t = t.sort_values(["node_id", "start_day"])
t.to_csv(os.path.join(d, "ITN_v3.csv"), index=False)
print("rows", len(t))
print(t[["node_id", "year", "start_day", "coverage"]].to_string(index=False))
