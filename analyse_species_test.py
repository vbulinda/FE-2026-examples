import glob
import pandas as pd
ids = {1.0: "22f2e1b8-d5ee-4034-9cdc-d1f12ca87c22", 0.5: "c5491159-91a9-4e75-9576-3a51a47cbf54",
       0.25: "fa5bbb57-d2bb-4a6f-91ff-5c749d7843de", 0.125: "f8903f32-0c63-4d1c-8192-208ced019329",
       0.0625: "0257defb-3920-40de-a239-f96663bba8b0"}
for m, i in ids.items():
    f = glob.glob(f"experiments/*{i}*/*/output/ReportVectorStats.csv")
    if not f:
        print(m, "file not found")
        continue
    d = pd.read_csv(f[0])
    if m == 1.0:
        print("columns:", list(d.columns))
    if "VectorPopulation" not in d.columns:
        print(m, "no VectorPopulation column")
        continue
    d = d[d["Time"] >= 730]
    s = d.groupby("Species")["VectorPopulation"].sum()
    p = (100 * s / s.sum()).round(1)
    print("ARA_MULT", m, ":", ", ".join(k + " " + str(v) + "%" for k, v in p.items()))
