import os, datetime as dt, pandas as pd
d = os.path.join("inputs", "interventions")
t = pd.read_csv(os.path.join(d, "ITN_v4.csv"))
day = (dt.date(2015, 12, 1) - dt.date(2010, 1, 1)).days
for n in [5, 6, 7]:
    r = t[(t.node_id == n) & (t.year.isin([2014, 2015]))]
    if r.empty:
        print("node", n, "no 2014/2015 row"); continue
    k = r.sort_values("year").index[-1]
    t = t.drop(index=r.index.difference([k]))
    t.loc[k, "start_day"] = day
    t.loc[k, "year"] = 2015
t = t.sort_values(["node_id", "start_day"])
t.to_csv(os.path.join(d, "ITN_v5.csv"), index=False)
print("start_day for 1 Dec 2015 =", day)
print(t[["node_id", "year", "start_day", "coverage"]].to_string(index=False))
