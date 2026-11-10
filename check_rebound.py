import glob, json
import numpy as np, pandas as pd
ids = {"imp0": "ba73cdf6-990e-4b85-9248-006a178dbe10", "imp1": "8080fc65-416e-4a79-8260-68e2cf84d55f",
       "imp3": "e81498fb-4146-4318-bb52-7b7738270df1", "imp10": "4e4c08c8-8a44-4a64-8d6a-dcaacf80b3b0"}
rows = []
for sc, i in ids.items():
    for cfg in glob.glob(f"experiments/*{i}*/*/config.json"):
        c = json.load(open(cfg))["parameters"]
        m = json.load(open(cfg.replace("config.json", "output/MalariaSummaryReport_monthly.json")))
        t = np.array(m["DataByTime"]["Time Of Report"], dtype=float)
        pf = np.array(m["DataByTimeAndAgeBins"]["PfPR by Age Bin"], dtype=float)[:, 1] * 100
        for a, b in zip(t, pf):
            if a >= 3400:
                rows.append((sc, round(c["x_Temporary_Larval_Habitat"], 2), int(a), b))
d = pd.DataFrame(rows, columns=["scenario", "xTLH", "day", "pfpr"])
for x in sorted(d["xTLH"].unique()):
    tab = d[d["xTLH"] == x].groupby(["day", "scenario"])["pfpr"].median().unstack()[list(ids)].round(1)
    print(f"\nxTLH {x}: PfPR children 0.25-5 y (%), median of 2 seeds, by report day (survey window = days 3961-4006)")
    print(tab.to_string())
