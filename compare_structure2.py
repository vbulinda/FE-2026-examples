import glob, json
import numpy as np, pandas as pd
ids = {"d150_imp1": "94b4f782-b73f-4222-9b24-82449d9afd9a", "d150_imp3": "80507db9-666a-41f4-a873-305c6306652a",
       "d150_imp10": "40ed82f7-d52d-4e4a-b730-02cd4a4c77e7", "d300_imp1": "5ba92d30-c2b0-4642-b560-7d31a3eebd61",
       "d300_imp3": "bccfa802-781a-4ccf-b24a-c0a0911b35bf", "d300_imp10": "9dc4dd6b-88cd-485b-82f7-0fc1925243e8"}
order = list(ids)
rows = []
for sc, i in ids.items():
    for cfg in glob.glob(f"experiments/*{i}*/*/config.json"):
        c = json.load(open(cfg))["parameters"]
        m = json.load(open(cfg.replace("config.json", "output/MalariaSummaryReport_monthly.json")))
        t = np.array(m["DataByTime"]["Time Of Report"], dtype=float)
        pf = np.array(m["DataByTimeAndAgeBins"]["PfPR by Age Bin"], dtype=float)[:, 1] * 100
        for a, b in zip(t, pf):
            rows.append((sc, round(c["x_Temporary_Larval_Habitat"], 2), c["Run_Number"], int(a), b))
d = pd.DataFrame(rows, columns=["scenario", "xTLH", "seed", "day", "pfpr"])
print("simulations per scenario:", d.groupby("scenario")["seed"].nunique().to_dict())
w = d[(d["day"] >= 3945) & (d["day"] <= 4020)].groupby(["scenario", "xTLH"])["pfpr"].median().unstack().loc[order].round(1)
print("\nPfPR children 0.25-5 y (%) in the survey window (reports of days 3961-4015), median of seeds and reports:")
print(w.to_string())
for x in sorted(d["xTLH"].unique()):
    tab = d[(d["xTLH"] == x) & (d["day"] >= 3600)].groupby(["day", "scenario"])["pfpr"].median().unstack()[order].round(1)
    print(f"\nxTLH {x}: monthly PfPR (%), median of 2 seeds")
    print(tab.to_string())
print("\nSurvey, Migori children 6-59 months, Nov-Dec 2020: microscopy 13.6% (4.7-21.9), RDT 17.4% (7.5-28.1)")
