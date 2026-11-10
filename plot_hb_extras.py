"""
plot_hb_extras.py - Homa Bay (node 4) seasonal profile and clinical incidence by year (model output, NOT calibrated).
    python plot_hb_extras.py <output_folder_name> [xTLH_main] [xTLH_other1 xTLH_other2]
e.g. python plot_hb_extras.py hb_A_s10_k5 4.0 3.1748 5.0397
Reads experiments/my_outputs/<folder>/SpatialReportMalariaFiltered_children.csv and _all_ages.csv, and the NASA POWER rainfall file
inputs/climate/LAKE_HOMABAY_NASA/2010001-2026304/dtk_15arcmin_rainfall_daily.bin (float32, one value per day, day 0 = 1 Jan 2010, 365-day years).
Writes homabay_seasonal.png, homabay_clinical.png and the numbers behind them (homabay_seasonal.csv, homabay_clinical_by_year.csv) in the current folder.
Seasonal profile uses 2010-2018 only (before the 2019-2020 IRS rounds). Model prevalence = raw channel / 0.1635 (as in the calibration).
Clinical incidence: New_Clinical_Cases per year per 1000 people, whole population (the filtered reports do not age-filter this channel); 2026 is partial and left out.
There are no clinical-incidence data in the calibration: these panels are model output only.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FOLDER = sys.argv[1]
XMAIN = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
XOTH = [float(v) for v in sys.argv[3:5]] if len(sys.argv) > 3 else []
NODE, CONST = 4, 0.1635
OUT = os.path.join("experiments", "my_outputs", FOLDER)
RAIN = "inputs/climate/LAKE_HOMABAY_NASA/2010001-2026304/dtk_15arcmin_rainfall_daily.bin"
INK, MUTE, GRID, RED = "#1f2933", "#6b7785", "#e3e7eb", "#d7191c"
LIGHT = ["#ef8a8a", "#7a0a0c"]
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
MLEN = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MSTART = np.concatenate([[0], np.cumsum(MLEN)])
plt.rcParams.update({"font.size": 10, "text.color": INK, "axes.labelcolor": INK, "xtick.color": MUTE, "ytick.color": MUTE})


def clean(ax):
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(MUTE)


def load(name):
    d = pd.read_csv(os.path.join(OUT, f"SpatialReportMalariaFiltered_{name}.csv"))
    d = d[d["node"] == NODE].copy()
    d["xTLH"] = d["xTLH"].astype(float).round(4)
    d["year"] = 2010 + (d["time"] // 365).astype(int)
    d["doy"] = (d["time"] % 365).astype(int)
    return d


ch, al = load("children"), load("all_ages")
for nm, d in (("children", ch), ("all_ages", al)):
    dt = np.diff(np.sort(d["time"].unique()))
    print(f"{nm}: {d['Run_Number'].nunique()} seeds, xTLH {sorted(d['xTLH'].unique())}, time {d['time'].min()}-{d['time'].max()}, step {np.median(dt):.0f} d (min {dt.min():.0f}, max {dt.max():.0f})")
xs_all = [XMAIN] + XOTH
for x in xs_all:
    if x not in set(ch["xTLH"].round(4)):
        raise SystemExit(f"xTLH {x} not in the output; available {sorted(ch['xTLH'].unique())}")

# ---- seasonal profile, 2010-2018 ----
rain = np.fromfile(RAIN, dtype="float32")
ry = rain[: 9 * 365].reshape(9, 365)
rain_month = np.array([ry[:, MSTART[m]:MSTART[m + 1]].sum(axis=1).mean() for m in range(12)])   # mean mm per month over 2010-2018
pre = ch[ch["year"] <= 2018].copy()
pre["month"] = np.searchsorted(MSTART[1:], pre["doy"], side="right")
pre["p"] = pre["PCR_Parasite_Prevalence"] / CONST * 100
rows = []
for x in xs_all:
    g = pre[pre["xTLH"] == x].groupby(["Run_Number", "month"])["p"].mean().unstack()   # seed x month
    for m in range(12):
        rows.append(dict(xTLH=x, month=MONTHS[m], median=g[m].median(), min=g[m].min(), max=g[m].max(), rain_mm=rain_month[m]))
S = pd.DataFrame(rows); S.to_csv("homabay_seasonal.csv", index=False)
fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7.2), sharex=True, gridspec_kw={"height_ratios": [1, 1.5]})
clean(a1); a1.bar(range(12), rain_month, color="#a9c7dc", width=0.7)
a1.set_ylabel("Rainfall (mm per month)"); a1.set_title("A. Rainfall at Homa Bay, mean of 2010-2018 (NASA POWER)", loc="left", fontsize=11, fontweight="bold")
clean(a2)
for i, x in enumerate(XOTH):
    s = S[S["xTLH"] == x]; a2.plot(range(12), s["median"], color=LIGHT[i], lw=1.4, label=f"xTLH {x:.2f}")
s = S[S["xTLH"] == XMAIN]
a2.fill_between(range(12), s["min"], s["max"], color=RED, alpha=0.15, lw=0)
a2.plot(range(12), s["median"], color=RED, lw=2.4, marker="o", ms=5, label=f"xTLH {XMAIN:.2f} (median of seeds; band = seed range)")
a2.set_ylabel("Model PCR prevalence, children 6-59 months (%)"); a2.set_ylim(bottom=0)
a2.set_xticks(range(12)); a2.set_xticklabels(MONTHS)
a2.set_title("B. Model prevalence by calendar month, 2010-2018 (before IRS)", loc="left", fontsize=11, fontweight="bold")
a2.legend(frameon=False, fontsize=8.5, loc="lower left")
fig.suptitle("Homa Bay County (node 4): seasonal profile, model output", fontsize=13, fontweight="bold", x=0.01, ha="left")
fig.tight_layout(rect=(0, 0.03, 1, 0.96))
fig.text(0.01, 0.005, "Seasonal timing is not calibrated: the fit uses two survey points (MIS 2015, 2020). 365-day model years, day 0 = 1 Jan.", fontsize=8.5, color=MUTE)
fig.savefig("homabay_seasonal.png", dpi=200, bbox_inches="tight")

# ---- clinical incidence per year ----
# New_Clinical_Cases and Population in the filtered reports are whole-node values (identical in the children and all-ages files),
# so this is population-wide incidence, NOT children only.
def seeds(d, x):
    d = d[(d["xTLH"] == x) & (d["year"] <= 2025)]
    g = d.groupby(["Run_Number", "year"]).agg(cases=("New_Clinical_Cases", "sum"), pop=("Population", "mean")).reset_index()
    g["inc"] = g["cases"] / g["pop"] * 1000
    return g.pivot(index="year", columns="Run_Number", values="inc")
R = []
fig, ax = plt.subplots(figsize=(9, 5))
clean(ax); ax.axvspan(2018.5, 2020.5, color="#e1edf4", zorder=0); ax.axvspan(2020.5, 2025.5, color="#f1f3f5", zorder=0)
for i, x in enumerate(XOTH):
    t = seeds(al, x).median(axis=1); ax.plot(t.index, t.values, color=LIGHT[i], lw=1.4, label=f"xTLH {x:.2f} (median)")
P = seeds(al, XMAIN)
for c in P.columns: ax.plot(P.index, P[c], color=RED, lw=0.7, alpha=0.35)
ax.plot(P.index, P.median(axis=1), color=RED, lw=2.6, marker="o", ms=5, label=f"xTLH {XMAIN:.2f} (median of seeds; thin lines = each seed)")
ax.set_ylim(bottom=0); ax.set_xlabel("Year"); ax.set_ylabel("Clinical cases per 1000 per year")
ax.legend(frameon=False, fontsize=8.5, loc="upper left")
for yr in P.index: R.append(dict(year=yr, median=P.loc[yr].median(), min=P.loc[yr].min(), max=P.loc[yr].max()))
fig.suptitle("Homa Bay County (node 4): clinical incidence by year, model output (not calibrated)", fontsize=12.5, fontweight="bold", x=0.01, ha="left")
fig.tight_layout(rect=(0, 0.07, 1, 0.95))
fig.text(0.01, 0.005, "Whole-population counts (the report does not separate children). No clinical data used. Blue band: IRS 2019-2020; grey: after 2020, further IRS not modelled,\nso the 2021 rebound is partly an artefact and seeds diverge. Case management: lake-zone values. 2026 omitted (partial year).", fontsize=8, color=MUTE)
fig.savefig("homabay_clinical.png", dpi=200, bbox_inches="tight")
C = pd.DataFrame(R); C.to_csv("homabay_clinical_by_year.csv", index=False)
print("\nSeasonal (xTLH %.2f): prevalence by month, median %%, and rainfall mm" % XMAIN)
print(S[S["xTLH"] == XMAIN][["month", "median", "min", "max", "rain_mm"]].round(1).to_string(index=False))
print("\nClinical cases per 1000 people per year, whole population (xTLH %.2f), median [min-max over seeds]" % XMAIN)
print("; ".join(f"{int(r.year)} {r['median']:.0f} [{r['min']:.0f}-{r['max']:.0f}]" for _, r in C.iterrows()))
print("\nSaved: homabay_seasonal.png homabay_clinical.png homabay_seasonal.csv homabay_clinical_by_year.csv")
