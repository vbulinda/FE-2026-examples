"""
score_bungoma.py - score the Bungoma pickup against Kenya MIS 2015 and 2020 (children 6-59 months).
    python score_bungoma.py [path_to_children_csv]
Changes from score_calibration14.py (see EMOD_Bungoma_Calibration_Log.md, step 17):
  - Bungoma (node 3) only
  - model prevalence taken over the survey FIELDWORK days, not the whole calendar year:
      MIS 2015: 6 Jul-15 Aug 2015 = model days 2011-2051;  MIS 2020: 8 Nov-23 Dec 2020 = model days 3961-4006
  - RDT survey is compared with model PCR; microscopy survey with model Blood_Smear (not PCR)
  - model prevalence = raw channel / 0.1635 (share of the node aged 6-59 months, measured from the model, log step 16c)
  - J is computed per seed, then summarised by the median over seeds (extinct and persisting seeds do not average out)
  - SE_model per seed = binomial SE for the children in the node (0.1635 x 5000)
J(seed) = sum over survey years [p_model - p_survey]^2 / (SE_survey^2 + SE_model^2)
"""
import sys
import numpy as np, pandas as pd

NODE, COUNTY = 7, "bungoma"
CONST = 0.1635
N_CHILD = CONST * 5000
WINDOWS = {2015: (2011, 2051), 2020: (3961, 4006)}
CHANNEL = {"RDT": "PCR_Parasite_Prevalence", "microscopy": "Blood_Smear_Parasite_Prevalence"}
cpath = sys.argv[1] if len(sys.argv) > 1 else "experiments/my_outputs/lake_bungoma_pickup_v2/SpatialReportMalariaFiltered_children.csv"

sp = pd.read_csv(cpath)
sp = sp[sp["node"] == NODE].copy()
sp["xTLH"] = sp["xTLH"].astype(float).round(4)
mis = pd.read_csv("mis_prevalence_6_59m_by_county.csv")
mis = mis[mis["county"].str.lower() == COUNTY].copy()
mis["p"] = mis["prev_pct"] / 100
mis["se"] = ((mis["ci_hi"] - mis["ci_lo"]) / (2 * 1.96) / 100).clip(lower=0.01)

rows = []
for test, ch in CHANNEL.items():
    for y, (lo, hi) in WINDOWS.items():
        w = sp[(sp["time"] >= lo) & (sp["time"] <= hi)]
        g = w.groupby(["xTLH", "Run_Number"])[ch].mean().reset_index()
        g["p_model"] = (g[ch] / CONST).clip(upper=1.0)
        s = mis[(mis["test"] == test) & (mis["year"] == y)].iloc[0]
        g["test"], g["year"], g["p_obs"], g["se_obs"] = test, y, s["p"], s["se"]
        pm = g["p_model"].clip(0.001, 0.999)
        g["se_model"] = np.sqrt(pm * (1 - pm) / N_CHILD)
        g["z2"] = (g["p_model"] - g["p_obs"]) ** 2 / (g["se_obs"] ** 2 + g["se_model"] ** 2)
        rows.append(g[["test", "year", "xTLH", "Run_Number", "p_model", "p_obs", "se_obs", "se_model", "z2"]])
d = pd.concat(rows)
d.to_csv("bungoma_scores_per_seed_year.csv", index=False)
J = d.groupby(["test", "xTLH", "Run_Number"])["z2"].sum().reset_index(name="J")
J.to_csv("bungoma_J_per_seed.csv", index=False)

print("Survey targets (children 6-59 months), Bungoma:")
print(mis[["year", "test", "prev_pct", "ci_lo", "ci_hi", "n_tested"]].to_string(index=False))
for test in CHANNEL:
    t = J[J["test"] == test]
    s = t.groupby("xTLH")["J"].agg(["median", "mean", "min", "max"]).round(1)
    p = (d[d["test"] == test].groupby(["xTLH", "year"])["p_model"].median().unstack() * 100).round(1)
    p.columns = [f"model{c}%" for c in p.columns]
    print(f"\n=== Bungoma, model {CHANNEL[test]} vs survey {test}: J by xTLH over the 4 seeds (model % = median over seeds) ===")
    print(s.join(p).to_string())
    best = s["median"].idxmin()
    ok = s[s["median"] <= s["median"].min() + 3.84].index
    print(f"Lowest median J: xTLH {best} (J={s['median'].min():.1f}); xTLH within 3.84 of it: {ok.min()} to {ok.max()}")
    if best in (s.index.min(), s.index.max()):
        print("WARNING: the lowest J is at the edge of the xTLH grid; the true best may lie outside it.")
print("\nSaved: bungoma_scores_per_seed_year.csv, bungoma_J_per_seed.csv")
print("Reminder: two survey points and one parameter; J near 0-4 = consistent. Survey CIs are wide (n=127 and 242).")
