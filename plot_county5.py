"""
plot_county.py - calibration figure for one county from the EMOD experiment outputs.
    python plot_county.py <county> <output_folder_name> [xTLH1 xTLH2 xTLH3]
    e.g.  python plot_county.py bungoma lake_bungoma_v2 1.3348 1.7818 3.1748
Reads (written by analyzer_migori2.py and score_<county>.py, run from FE-2026-examples):
    experiments/my_outputs/<folder>/SpatialReportMalariaFiltered_children.csv
    <county>_J_per_seed.csv, <county>_scores_per_seed_year.csv, mis_prevalence_6_59m_by_county.csv
Writes experiments/my_outputs/<folder>/<county>_calibration.png and .pdf
Model prevalence in children = raw spatial channel / 0.1635. Data in black, model in red.
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COUNTY, FOLDER = sys.argv[1].lower(), sys.argv[2]
SHOW = [float(x) for x in sys.argv[3:6]] or [1.3348, 1.7818, 3.1748]
OUT = os.path.join("experiments", "my_outputs", FOLDER)
CONST = 0.1635
TAG = os.environ.get("TAG", "")
LABEL = os.environ.get("LABEL", "")
SUF = f"_{TAG}" if TAG else ""
CHANNEL = {"RDT": "PCR_Parasite_Prevalence", "microscopy": "Blood_Smear_Parasite_Prevalence"}
EVENTS = {"bungoma": [(2160, "LLIN 2015"), (2800, "LLIN 2017"), (4180, "PBO 2021"), (5291, "PBO 2024")],
          }.get(COUNTY, [])
NETS = {"migori": [1730, 2661]}.get(COUNTY, [])   # LLIN campaigns 2014, 2017 (dotted lines)
IRSB = {"migori": [(2616, 300), (2981, 300), (3312, 300), (3690, 150)]}.get(COUNTY, [])   # IRS rounds: start day, full-effect days (sourced-rounds runs)
MIS = pd.read_csv("mis_prevalence_6_59m_by_county.csv")
MIS = MIS[MIS["county"].str.lower() == COUNTY]
NODE = int(MIS["node_id"].iloc[0])
YEARPOS = {2015: 2015.55, 2020: 2020.9}
SURVEY = {t: {YEARPOS[int(r.year)]: (r.prev_pct, r.ci_lo, r.ci_hi) for r in g.itertuples()}
          for t, g in MIS.groupby("test")}
INK, MUTE, GRID, RED = "#1f2933", "#6b7785", "#e3e7eb", "#d7191c"
SHADES = dict(zip(SHOW, ["#ef8a8a", "#d7191c", "#7a0a0c"]))
plt.rcParams.update({"font.size": 10, "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTE, "ytick.color": MUTE})

sp = pd.read_csv(os.path.join(OUT, "SpatialReportMalariaFiltered_children.csv"))
sp = sp[sp["node"] == NODE].copy()
sp["xTLH"] = sp["xTLH"].astype(float).round(4)
J = pd.read_csv(f"{COUNTY}_J_per_seed.csv" if not TAG else f"{COUNTY}_J_{TAG}.csv"); J["xTLH"] = J["xTLH"].astype(float).round(4)
Y = pd.read_csv(f"{COUNTY}_scores_per_seed_year.csv" if not TAG else f"{COUNTY}_scores_{TAG}.csv")
Y["xTLH"] = Y["xTLH"].astype(float).round(4)
SHOW = [min(sorted(J["xTLH"].unique()), key=lambda v: abs(v - x)) for x in SHOW]
SHADES = dict(zip(SHOW, ["#ef8a8a", "#d7191c", "#7a0a0c"]))


def clean(ax):
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(MUTE)


def series_panel(ax, test, title):
    clean(ax)
    top = 0
    ax.axvspan(2021, 2026.85, color="#f1f3f5", zorder=0)
    for x in SHOW:
        t = sp[sp["xTLH"] == x]
        w = t.pivot_table(index="time", columns="Run_Number", values=CHANNEL[test]) / CONST * 100
        w = w.rolling(30, center=True, min_periods=10).mean()
        yr = 2010 + (w.index.values - 1) / 365
        ax.fill_between(yr, w.min(axis=1), w.max(axis=1), color=SHADES[x], alpha=0.15, lw=0)
        ax.plot(yr, w.median(axis=1), color=SHADES[x], lw=1.6, label=f"xTLH {x:.2f}")
        top = max(top, float(w.max().max()))
    for yr, (p, lo, hi) in SURVEY[test].items():
        ax.errorbar([yr], [p], yerr=[[p - lo], [hi - p]], fmt="o", color="black",
                    ms=8, capsize=4, lw=1.8, zorder=5)
        top = max(top, hi)
    ax.plot([], [], "o", color="black", label=f"{test} survey (95% CI)")
    top = 10 * np.ceil((top + 20) / 10)
    for d, lab in EVENTS:
        ax.axvline(2010 + d / 365, color=MUTE, ls=":", lw=1)
        ax.text(2010 + d / 365 + 0.05, top * 0.97, lab, rotation=90,
                va="top", fontsize=8, color=MUTE)
    for i, (d, n) in enumerate(IRSB):
        ax.axvspan(2010 + d / 365, 2010 + (d + n) / 365, color="#9ecae1", alpha=0.3, lw=0, zorder=1,
                   label="IRS full-effect period" if i == 0 else None)
    for i, d in enumerate(NETS):
        ax.axvline(2010 + d / 365, color=MUTE, ls=":", lw=1.3, label="LLIN campaign" if i == 0 else None)
    ax.text(2026.9, top * 0.97, "not fitted (post-2020\ninterventions incomplete)", fontsize=8.5, color=MUTE, ha="right",
            va="top", zorder=6, bbox=dict(fc="none", ec="none"))
    ax.set_xlim(2010, 2027); ax.set_ylim(0, top)
    ax.set_xlabel("Year"); ax.set_ylabel("Prevalence, children 6-59 months (%)")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=True, facecolor="white", framealpha=0.92, edgecolor="none", fontsize=8.5, loc="upper left", ncol=2 if (IRSB or NETS) else 1)


fig, axs = plt.subplots(2, 2, figsize=(13, 9), dpi=200)
series_panel(axs[0, 0], "RDT", "A. Model PCR through time vs RDT survey")
series_panel(axs[0, 1], "microscopy", "B. Model blood-smear through time vs microscopy survey")

ax = axs[1, 0]; clean(ax)
for test, ls, mk, c in (("RDT", "-", "o", "black"), ("microscopy", "--", "s", RED)):
    t = J[J["test"] == test].groupby("xTLH")["J"].agg(["median", "min", "max"])
    ax.fill_between(t.index, t["min"], t["max"], color=c, alpha=0.12, lw=0)
    ax.plot(t.index, t["median"], color=c, ls=ls, marker=mk, ms=5, lw=2, label=f"vs {test}")
    b = t["median"].idxmin()
    ax.plot([b], [t.loc[b, "median"]], "o", ms=12, mfc="none", mec=INK, mew=1.4)
ax.axhline(4, color="#9aa5b1", lw=1, ls=":", label="J = 4: close fit")
lo_, hi_ = J["xTLH"].min() * 0.8, J["xTLH"].max() * 1.2
xt = [v for v in (0.5, 1, 2, 4, 8, 16, 32) if lo_ <= v <= hi_]
ax.set_xscale("log"); ax.set_xticks(xt)
ax.set_xticklabels([str(v) for v in xt]); ax.minorticks_off()
ax.set_ylim(-1.5, 30)
ax.set_xlabel("x_Temporary_Larval_Habitat (xTLH)")
ax.set_ylabel("J (lower = better fit)")
ax.set_title("C. Fit to both surveys (rings = lowest median J)",
             loc="left", fontsize=11, fontweight="bold")
ax.legend(frameon=False, fontsize=8.5, loc="lower right")

ax = axs[1, 1]; clean(ax)
r = Y[Y["test"] == "RDT"]
for i, x in enumerate(SHOW):
    for k, (yr, col) in enumerate(((2015, RED), (2020, "black"))):
        d = r[(r["xTLH"] == x) & (r["year"] == yr)]
        xs = i + (k - 0.5) * 0.3 + np.linspace(-0.05, 0.05, max(len(d), 1))[: len(d)]
        lab = f"model {yr}" if i == 0 else None
        ax.plot(xs, d["p_model"] * 100, "o", color=col, ms=6, alpha=0.8, label=lab)
for yr, col in ((2015, RED), (2020, "black")):
    p, lo, hi = SURVEY["RDT"][YEARPOS[yr]]
    ax.axhspan(lo, hi, color=col, alpha=0.07); ax.axhline(p, color=col, ls="--", lw=1.4, label=f"{yr} survey {p}%")
    ax.axhline(lo, color=col, ls=":", lw=0.8); ax.axhline(hi, color=col, ls=":", lw=0.8)
ax.set_xticks(range(len(SHOW))); ax.set_xticklabels([f"xTLH {x:.2f}" for x in SHOW])
ax.set_xlim(-0.6, len(SHOW) - 0.4); ax.set_ylim(0, 68)
ax.set_ylabel("Model PCR prevalence in the survey window (%)")
ax.set_title("D. Each seed vs the RDT survey (shaded: 95% CI)",
             loc="left", fontsize=11, fontweight="bold")
ax.legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, columnspacing=1.2, handlelength=1.6)

fig.suptitle(f"{COUNTY.capitalize()} County (node {NODE}): EMOD vs Kenya MIS 2015, 2020",
             fontsize=13, fontweight="bold", x=0.01, ha="left")
fig.tight_layout(rect=(0, 0.02, 1, 0.96))
fig.text(0.01, 0.005, LABEL, fontsize=8.5, color=MUTE)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"{COUNTY}_calibration{SUF}.{ext}"), bbox_inches="tight")
print("Saved:", os.path.join(OUT, f"{COUNTY}_calibration{SUF}.png"))
