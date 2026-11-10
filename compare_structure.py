import glob, json
import numpy as np, pandas as pd
ids = {"imp0": "ba73cdf6-990e-4b85-9248-006a178dbe10", "imp1": "8080fc65-416e-4a79-8260-68e2cf84d55f",
       "imp3": "e81498fb-4146-4318-bb52-7b7738270df1", "imp10": "4e4c08c8-8a44-4a64-8d6a-dcaacf80b3b0"}
order = list(ids)
rows = []
for sc, i in ids.items():
    for cfg in glob.glob(f"experiments/*{i}*/*/config.json"):
        c = json.load(open(cfg))["parameters"]
        out = cfg.replace("config.json", "output/")
        m = json.load(open(out + "MalariaSummaryReport_monthly.json"))
        t = np.array(m["DataByTime"]["Time Of Report"], dtype=float)
        eir = np.array(m["DataByTime"]["Annual EIR"], dtype=float)
        pf = np.array(m["DataByTimeAndAgeBins"]["PfPR by Age Bin"], dtype=float)[:, 1] * 100
        yr = 2010 + ((t - 1) // 365).astype(int)
        for y in np.unique(yr):
            rows.append((sc, round(c["x_Temporary_Larval_Habitat"], 2), c["Run_Number"], int(y),
                         eir[yr == y].mean(), pf[yr == y].mean()))
d = pd.DataFrame(rows, columns=["scenario", "xTLH", "seed", "year", "eir", "pfpr"])
print("simulations per scenario:", d.groupby("scenario")["seed"].count().div(11).to_dict())
for x in sorted(d["xTLH"].unique()):
    s = d[d["xTLH"] == x]
    for col, label in (("pfpr", "PfPR in children 0.25-5 y (%)"), ("eir", "Annual EIR")):
        tab = s.groupby(["year", "scenario"])[col].median().unstack()[order].round(1)
        print(f"\nxTLH {x}: {label}, median of 2 seeds")
        print(tab.to_string())
print("\nSurvey, Migori children 6-59 months: 2015 microscopy 31.1% / RDT 35.6%;  2020 microscopy 13.6% / RDT 17.4%")
