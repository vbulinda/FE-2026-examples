import sys
import numpy as np, pandas as pd
SHARE = 0.1615
NAMES = {1: "kisumu", 2: "siaya", 3: "migori", 4: "homa bay", 5: "kakamega", 6: "busia", 7: "bungoma", 8: "vihiga"}
sp = pd.read_csv(sys.argv[1]); sp["xTLH"] = sp["xTLH"].astype(float).round(4)
sp["p"] = (sp["PCR_Parasite_Prevalence"] / SHARE).clip(upper=1.0)
mis = pd.read_csv(sys.argv[2] if len(sys.argv) > 2 else "mis_prevalence_6_59m_by_county.csv"); mis = mis[mis["county"] != "ZONE"].copy()
mis["se"] = ((mis["ci_hi"] - mis["ci_lo"]) / (2 * 1.96) / 100).clip(lower=0.01); mis["pm"] = mis["prev_pct"] / 100
rows = []
for y in (2015, 2020):
    lo, hi = (y - 2010) * 365, (y - 2009) * 365
    g = sp[(sp.time >= lo) & (sp.time < hi)].groupby(["xTLH", "node", "Run_Number"])["p"].mean().reset_index(); g["year"] = y; rows.append(g)
s = pd.concat(rows); s["county"] = s["node"].map(NAMES)
ext = s[s.year == 2020].assign(ext=lambda d: d.p < 0.005).groupby(["county", "xTLH"])["ext"].agg(["sum", "count"]).reset_index()
for test in ("RDT", "microscopy"):
    m = s.merge(mis[mis.test == test][["county", "year", "pm", "se"]], on=["county", "year"]); m["z2"] = (m.p - m.pm) ** 2 / m.se ** 2
    Js = m.groupby(["county", "xTLH", "Run_Number"])["z2"].sum().reset_index(name="J")
    J = Js.groupby(["county", "xTLH"])["J"].median().reset_index()
    print(f"\n=== median-over-seeds J, MIS {test} ==="); print(J.pivot(index="county", columns="xTLH", values="J").round(1).to_string())
    print(f"\nBest xTLH per county ({test}):")
    tot = 0
    for c, gg in J.groupby("county"):
        b = gg.loc[gg.J.idxmin()]; e = ext[(ext.county == c) & (ext.xTLH == b.xTLH)].iloc[0]; tot += b.J
        det = m[(m.county == c) & (m.xTLH == b.xTLH)].groupby("year").agg(p=("p", "median"), pm=("pm", "first"))
        pr = "  ".join(f"{y}: median seed {r.p*100:.1f}% vs MIS {r.pm*100:.1f}%" for y, r in det.iterrows())
        print(f"  {c:9s} best {b.xTLH:.4f} (median J={b.J:.1f})  2020-extinct seeds {int(e['sum'])}/{int(e['count'])}   {pr}")
    print("sum of best median J:", round(tot, 1))
