import os, glob, json, argparse, datetime as dt
import numpy as np, pandas as pd
import manifest
START = dt.date(2010, 1, 1)
NAMES = {1: "kisumu", 2: "siaya", 3: "migori", 4: "homa bay", 5: "kakamega", 6: "busia", 7: "bungoma", 8: "vihiga"}
ap = argparse.ArgumentParser(); ap.add_argument("--exp", required=True); ap.add_argument("--name", required=True)
ap.add_argument("--mis", default="mis_prevalence_6_59m_by_county.csv"); ap.add_argument("--win2020", default="")
a = ap.parse_args()
def win(y, s):
    (m1, d1), (m2, d2) = [tuple(int(v) for v in p.split("-")) for p in s.split(":")]
    return (dt.date(y, m1, d1) - START).days, (dt.date(y, m2, d2) - START).days + 1
WIN = {"annual": {2015: win(2015, "01-01:12-31"), 2020: win(2020, "01-01:12-31")}}
WIN["months"] = {2015: win(2015, "07-06:08-15"), 2020: win(2020, a.win2020) if a.win2020 else WIN["annual"][2020]}
wdir = os.path.join(manifest.job_directory, "my_outputs", a.name)
cpath = (glob.glob(os.path.join(wdir, "**", "SpatialReportMalariaFiltered_children.csv"), recursive=True) or [None])[0]
if not cpath: raise SystemExit("children CSV not found: run the analyzer first")
shares = []
for f in glob.glob(os.path.join(manifest.job_directory, f"*{a.exp}*", "*", "output", "MalariaSummaryReport_monthly.json")):
    try:
        pop = np.array(json.load(open(f))["DataByTimeAndAgeBins"]["Average Population by Age Bin"])
        shares.append(((pop[:, 1] - pop[:, 0]) / pop.sum(axis=1)).mean())
    except Exception: pass
SHARE = float(np.mean(shares)) if shares else 0.1635
print(f"children share {SHARE:.4f} ({len(shares)} summary reports)" if shares else "WARNING: using share 0.1635")
sp = pd.read_csv(cpath); sp["xTLH"] = sp["xTLH"].astype(float).round(4)
sp["p"] = (sp["PCR_Parasite_Prevalence"] / SHARE).clip(upper=1.0)
mis = pd.read_csv(a.mis); mis = mis[mis["county"] != "ZONE"].copy()
mis["se"] = ((mis["ci_hi"] - mis["ci_lo"]) / (2 * 1.96) / 100).clip(lower=0.01); mis["pm"] = mis["prev_pct"] / 100
out = {}
for mode, W in WIN.items():
    rows = []
    for y, (lo, hi) in W.items():
        w = sp[(sp["time"] >= lo) & (sp["time"] < hi)]
        print(f"{mode} {y}: days {lo}-{hi}, {w['time'].nunique()} time points")
        g = w.groupby(["xTLH", "node", "Run_Number"])["p"].mean().reset_index(); g["year"] = y; rows.append(g)
    sim = pd.concat(rows); sim["county"] = sim["node"].map(NAMES)
    agg = sim.groupby(["xTLH", "county", "year"])["p"].agg(["mean", "std", "count"]).reset_index()
    agg["se_model"] = (agg["std"].fillna(0) / np.sqrt(agg["count"]))
    if mode == "months": agg.to_csv("calibration_model_prevalence_months.csv", index=False)
    res = []
    for test in ("RDT", "microscopy"):
        t = mis[mis["test"] == test][["county", "year", "pm", "se"]]
        m = agg.merge(t, on=["county", "year"]); m["z2"] = (m["mean"] - m["pm"]) ** 2 / (m["se"] ** 2 + m["se_model"] ** 2)
        J = m.groupby(["county", "xTLH"])["z2"].sum().reset_index(name="J"); J["test"] = test; res.append(J)
        if mode == "months":
            print(f"\n=== survey-month J by county and xTLH, MIS {test} ==="); print(J.pivot(index="county", columns="xTLH", values="J").round(1).to_string())
            print(f"\nBest xTLH per county ({test}); range = J within 3.84 of the minimum:")
            for c, gg in J.groupby("county"):
                b = gg.loc[gg["J"].idxmin()]; ok = gg[gg["J"] <= b["J"] + 3.84]["xTLH"]
                det = m[(m.county == c) & (m.xTLH == b["xTLH"])].sort_values("year")
                pr = "  ".join(f"{int(r.year)}: model {r['mean']*100:.1f}% vs MIS {r.pm*100:.1f}%" for _, r in det.iterrows())
                print(f"  {c:9s} best {b.xTLH:.4f} (J={b.J:.1f})  range {ok.min():.4f}-{ok.max():.4f}   {pr}")
    out[mode] = pd.concat(res)
    if mode == "months": out[mode].to_csv("calibration_J_by_county_xTLH_months.csv", index=False)
print("\n=== best J per county: annual mean vs survey months ===")
for test in ("RDT", "microscopy"):
    ba = out["annual"][out["annual"].test == test].groupby("county")["J"].min()
    bm = out["months"][out["months"].test == test].groupby("county")["J"].min()
    d = pd.DataFrame({"annual": ba, "survey_months": bm}).round(1); print(f"\n{test}"); print(d.to_string()); print("sum", d.sum().round(1).to_dict())
