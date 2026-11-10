import glob, json
import numpy as np, pandas as pd
EXP = "fbac22f7-9a53-45a1-8a2c-eff20584b415"
rows = []
for cfg in glob.glob(f"experiments/*{EXP}*/*/config.json"):
    c = json.load(open(cfg))["parameters"]
    out = cfg.replace("config.json", "output/")
    v = pd.read_csv(out + "ReportVectorStats.csv", usecols=["Time", "Species", "VectorPopulation"])
    v = v[(v["Time"] >= 1) & (v["Time"] <= 365)]
    s = v.groupby("Species")["VectorPopulation"].sum()
    m = json.load(open(out + "MalariaSummaryReport_monthly.json"))["DataByTime"]
    t = np.array(m["Time Of Report"], dtype=float); e = np.array(m["Annual EIR"], dtype=float)
    rows.append((round(c["x_Temporary_Latval_Habitat" if False else "x_Temporary_Larval_Habitat"], 2), c["Run_Number"],
                 100 * s.get("funestus", 0) / s.sum(), 100 * s.get("arabiensis", 0) / s.sum(), 100 * s.get("gambiae", 0) / s.sum(),
                 e[t <= 365].mean(), e[(t > 1825) & (t <= 2190)].mean()))
d = pd.DataFrame(rows, columns=["xTLH", "seed", "funestus%", "arabiensis%", "gambiae%", "EIR_2010", "EIR_2015"])
print(d.groupby("xTLH")[["funestus%", "arabiensis%", "gambiae%", "EIR_2010", "EIR_2015"]].median().round(1).to_string())
