import pandas as pd
d = pd.read_csv("experiments/my_outputs/lake_migori_pickup_v1/SpatialReportMalariaFiltered_children.csv")
d = d[d["node"] == 3].copy()
d["xTLH"] = d["xTLH"].astype(float).round(4)
d["year"] = 2010 + (d["time"] - 1) // 365
for x in (5.0397, 8.9797, 16.0):
    t = d[d["xTLH"] == x]
    g = t.groupby(["year", "Run_Number"])["PCR_Parasite_Prevalence"].mean().reset_index()
    g["pct"] = g["PCR_Parasite_Prevalence"] / 0.1635 * 100
    tab = g.pivot(index="year", columns="Run_Number", values="pct").round(1)
    tab["median"] = tab.median(axis=1).round(1)
    print(f"\nxTLH {x}: yearly mean PCR prevalence in children (%), one column per seed")
    print(tab.to_string())
