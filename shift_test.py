"""
shift_test.py - sensitivity of the Migori/county fit to the timing of the scoring windows (season timing).
    python shift_test.py <children_csv> <county> <node_id>
Both survey windows (MIS 2015: days 2011-2051, MIS 2020: days 3961-4006) are moved by -42..+42 days, J is
recomputed per seed and the median over seeds is reported per xTLH. Same J as score_county.py. Needs no new runs.
"""
import sys
import numpy as np, pandas as pd
cpath, county, node = sys.argv[1], sys.argv[2].lower(), int(sys.argv[3])
C = 0.1635; N = C * 5000
CH = {"RDT": "PCR_Parasite_Prevalence", "microscopy": "Blood_Smear_Parasite_Prevalence"}
W = {2015: (2011, 2051), 2020: (3961, 4006)}
sp = pd.read_csv(cpath); sp = sp[sp["node"] == node].copy(); sp["xTLH"] = sp["xTLH"].astype(float).round(4)
m = pd.read_csv("mis_prevalence_6_59m_by_county.csv"); m = m[m["county"].str.lower() == county].copy()
m["p"] = m["prev_pct"] / 100; m["se"] = ((m["ci_hi"] - m["ci_lo"]) / (2 * 1.96) / 100).clip(lower=0.01)
print("shift_days test best_xTLH minJ range_lo range_hi model2020%_at_best")
for sh in (-42, -28, -14, 0, 14, 28, 42):
    for test, ch in CH.items():
        parts = []
        for y, (lo, hi) in W.items():
            w = sp[(sp["time"] >= lo + sh) & (sp["time"] <= hi + sh)]
            g = w.groupby(["xTLH", "Run_Number"])[ch].mean().reset_index()
            g["pm"] = (g[ch] / C).clip(upper=1.0)
            s = m[(m["test"] == test) & (m["year"] == y)].iloc[0]
            q = g["pm"].clip(0.001, 0.999)
            g["z2"] = (g["pm"] - s["p"]) ** 2 / (s["se"] ** 2 + q * (1 - q) / N)
            g["year"] = y
            parts.append(g)
        d = pd.concat(parts)
        J = d.groupby(["xTLH", "Run_Number"])["z2"].sum().groupby("xTLH").median()
        ok = J[J <= J.min() + 3.84].index
        p20 = d[(d["year"] == 2020)].groupby("xTLH")["pm"].median()[J.idxmin()] * 100
        print(sh, test, J.idxmin(), round(J.min(), 1), ok.min(), ok.max(), round(p20, 1))
