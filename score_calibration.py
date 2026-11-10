"""
score_calibration.py - score the calibration pickup against Kenya MIS 2015 and 2020 (children 6-59 months).

    python score_calibration.py [path_to_children_csv] [path_to_mis_prevalence_csv]

Defaults: the children CSV written by analyzer_pickup3.py and mis_prevalence_6_59m_by_county.csv in the current folder.

For each county c and xTLH value th:   J(c,th) = sum over survey years y of [ p_model(c,y,th) - p_MIS(c,y) ]^2 / (SE_MIS^2 + SE_model^2)
p_model = mean PCR prevalence in children 6-59 months over the survey calendar year (pickup starts 1 Jan 2010), averaged over seeds.
SE_MIS from the survey 95% CI; SE_model from the seed spread (sd / sqrt(n seeds)); SE floor 0.01.
It scores against both RDT and microscopy because model PCR sits closest to RDT. Both results are printed; do not pick the better one afterwards without saying so.
"""
import sys, os, glob
import numpy as np, pandas as pd
import manifest
START_YEAR = 2010
NAMES = {1: "kisumu", 2: "siaya", 3: "migori", 4: "homa bay", 5: "kakamega", 6: "busia", 7: "bungoma", 8: "vihiga"}
SURVEYS = (2015, 2020)

wdir = os.path.join(manifest.job_directory, "my_outputs", "lake_pickup_cal")
cpath = sys.argv[1] if len(sys.argv) > 1 else (glob.glob(os.path.join(wdir, "**", "SpatialReportMalariaFiltered_children.csv"), recursive=True) or [None])[0]
mpath = sys.argv[2] if len(sys.argv) > 2 else "mis_prevalence_6_59m_by_county.csv"
if not cpath or not os.path.exists(cpath): raise SystemExit("children CSV not found: run analyzer_pickup3.py first")
sp = pd.read_csv(cpath); sp["xTLH"] = sp["xTLH"].astype(float).round(4)
mis = pd.read_csv(mpath); mis = mis[mis["county"] != "ZONE"].copy()
mis["se"] = ((mis["ci_hi"] - mis["ci_lo"]) / (2 * 1.96) / 100).clip(lower=0.01)
mis["p"] = mis["prev_pct"] / 100

# model prevalence per simulation, county and survey year
rows = []
for y in SURVEYS:
    lo, hi = (y - START_YEAR) * 365, (y - START_YEAR + 1) * 365
    w = sp[(sp["time"] >= lo) & (sp["time"] < hi)]
    g = w.groupby(["xTLH", "node", "Run_Number"])["PCR_Parasite_Prevalence"].mean().reset_index()
    g["year"] = y; rows.append(g)
sim = pd.concat(rows)
sim["county"] = sim["node"].map(NAMES)
agg = sim.groupby(["xTLH", "county", "year"])["PCR_Parasite_Prevalence"].agg(["mean", "std", "count"]).reset_index()
agg["se_model"] = (agg["std"].fillna(0) / np.sqrt(agg["count"])).clip(lower=0)
agg.to_csv("calibration_model_prevalence.csv", index=False)

res = []
for test in ("RDT", "microscopy"):
    t = mis[mis["test"] == test][["county", "year", "p", "se"]]
    m = agg.merge(t, on=["county", "year"])
    m["z2"] = (m["mean"] - m["p"]) ** 2 / (m["se"] ** 2 + m["se_model"] ** 2)
    J = m.groupby(["county", "xTLH"])["z2"].sum().reset_index(name="J"); J["test"] = test
    res.append(J)
    piv = J.pivot(index="county", columns="xTLH", values="J").round(1)
    print(f"\n=== J by county and xTLH, scored against MIS {test} ===")
    print(piv.to_string())
    print(f"\nBest xTLH per county ({test}); range = xTLH with J within 3.84 of the minimum (rough 95% band, grid only):")
    for c, gg in J.groupby("county"):
        b = gg.loc[gg["J"].idxmin()]; ok = gg[gg["J"] <= b["J"] + 3.84]["xTLH"]
        det = m[(m.county == c) & (m.xTLH == b["xTLH"])].sort_values("year")
        pr = "  ".join(f"{int(r.year)}: model {r['mean']*100:.1f}% vs MIS {r.p*100:.1f}%" for _, r in det.iterrows())
        print(f"  {c:9s} best {b.xTLH:.4f} (J={b.J:.1f}, 2 points)  range {ok.min():.4f}-{ok.max():.4f}   {pr}")
pd.concat(res).to_csv("calibration_J_by_county_xTLH.csv", index=False)
print("\nSaved: calibration_model_prevalence.csv, calibration_J_by_county_xTLH.csv")
print("Reminder: two survey points per county and one parameter; J near 0-4 means consistent, large J means the model cannot reproduce both surveys with this parameter alone.")