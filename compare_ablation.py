import glob, json
import numpy as np, pandas as pd
ids = {"full": "4f71ed74-62ec-4525-bdb9-27deb2e4c532", "no_irs": "2189f92d-d632-4f69-af83-8513f36ca08c",
       "no_itn17": "188591fd-49a6-4042-a07f-fa9d14da4456", "no_irs_itn17": "1c87bd32-f164-451e-af0c-0c241cae9087"}
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
