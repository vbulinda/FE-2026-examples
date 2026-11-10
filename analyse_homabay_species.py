"""
analyse_homabay_species.py - species mix of the 27 Homa Bay species-grid runs (year 3, Time >= 730, xTLH 4.0).
Run from FE-2026-examples after all 27 jobs finish:  python analyse_homabay_species.py
Reads homabay_species_ids.txt (lines g,a,f,expid) and each ReportVectorStats.csv; shares are of the summed VectorPopulation.
Mix A target: funestus about 85 percent. Mix B target: gambiae 32, arabiensis 49, funestus 19 percent.
"""
import glob
import pandas as pd

SPECIES = ["gambiae", "arabiensis", "funestus"]
rows = []
for line in open("homabay_species_ids.txt"):
    if not line.strip():
        continue
    g, a, f, eid = line.strip().split(",")
    fs = glob.glob(f"experiments/*{eid}*/**/ReportVectorStats.csv", recursive=True)
    if not fs:
        print(f"g {g} a {a} f {f}: not finished or file not found"); continue
    df = pd.read_csv(fs[0]); df = df[df["Time"] >= 730]
    by = df.groupby("Species")["VectorPopulation"].sum()
    pop = {s: float(by.get(s, 0.0)) for s in SPECIES}
    tot = sum(pop.values())
    if tot <= 0:
        print(f"g {g} a {a} f {f}: zero vectors"); continue
    sh = {s: 100 * pop[s] / tot for s in SPECIES}
    rows.append(dict(g=float(g), a=float(a), f=float(f), gam=sh["gambiae"], ara=sh["arabiensis"], fun=sh["funestus"],
                     vec_per_day=tot / df["Time"].nunique(), expid=eid))
r = pd.DataFrame(rows)
if r.empty:
    raise SystemExit("no finished runs")
r["distA"] = (r["fun"] - 85).abs()
r["distB"] = (r["gam"] - 32).abs() + (r["ara"] - 49).abs() + (r["fun"] - 19).abs()
pd.set_option("display.width", 200)
print(f"{len(r)} of 27 runs read")
print(r.drop(columns=["expid", "distA", "distB"]).round(1).to_string(index=False))
print("\nClosest to Mix A (funestus about 85):")
print(r.sort_values("distA").head(3)[["g", "a", "f", "gam", "ara", "fun"]].round(1).to_string(index=False))
print("\nClosest to Mix B (gam 32, ara 49, fun 19):")
print(r.sort_values("distB").head(3)[["g", "a", "f", "gam", "ara", "fun"]].round(1).to_string(index=False))
