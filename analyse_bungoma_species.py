"""
analyse_bungoma_species.py - species mix of the five Bungoma species-test runs (year 3, Time >= 730, xTLH 4.0).
Run from FE-2026-examples:  python analyse_bungoma_species.py
Reads ReportVectorStats.csv of each experiment; shares are of the summed VectorPopulation by species.
"""
import glob, sys
import pandas as pd

IDS = {5: "ca24a45b-c552-4d4c-8909-325ece4a28c3", 10: "1380ca6d-0186-4a5d-a816-c81a1410b9d2",
       20: "d6ef2284-6747-4cf0-a889-1c692a7d525c", 30: "840bdcd1-181f-486e-a184-02f4da54827c",
       40: "19f23a8d-7588-4362-97f4-99c78934ae77"}
SPECIES = ["gambiae", "arabiensis", "funestus"]
print("G_MULT |  gambiae % | arabiensis % | funestus % | mean vectors/day (gam, ara, fun)")
for g, eid in IDS.items():
    fs = glob.glob(f"experiments/*{eid}*/**/ReportVectorStats.csv", recursive=True)
    if not fs:
        print(f"{g:6d} | not finished or file not found"); continue
    df = pd.read_csv(fs[0]); df = df[df["Time"] >= 730]
    if "Species" in df.columns:
        by = df.groupby("Species")["VectorPopulation"].sum()
        pop = {s: float(by.get(s, 0.0)) for s in SPECIES}
        mean = {s: float(df[df["Species"] == s]["VectorPopulation"].mean()) for s in SPECIES}
    else:
        cols = {s: [c for c in df.columns if "VectorPopulation" in c and s in c] for s in SPECIES}
        if not all(cols.values()):
            print(f"{g:6d} | unknown layout; columns: {list(df.columns)}"); continue
        pop = {s: float(df[cols[s][0]].sum()) for s in SPECIES}
        mean = {s: float(df[cols[s][0]].mean()) for s in SPECIES}
    tot = sum(pop.values())
    print(f"{g:6d} | {100*pop['gambiae']/tot:11.1f} | {100*pop['arabiensis']/tot:12.1f} | {100*pop['funestus']/tot:10.1f} | "
          f"{mean['gambiae']:.0f}, {mean['arabiensis']:.0f}, {mean['funestus']:.0f}")
