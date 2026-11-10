"""
plot_migori.py - calibration figure for Migori from the EMOD experiment outputs.
    python plot_migori.py [tag]          tag = c067 (default), d150 or d225
Reads (all written by analyzer_migori2.py and score_migori.py):
    experiments/my_outputs/lake_migori_<tag>/SpatialReportMalariaFiltered_children.csv
    J_per_seed_<tag>.csv, scores_per_seed_year_<tag>.csv
Writes experiments/my_outputs/lake_migori_<tag>/Migori_calibration_<tag>.png and .pdf
Model prevalence in children = raw spatial channel / 0.1635 (share of the node aged 6-59 months).
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAG = sys.argv[1] if len(sys.argv) > 1 else "c067"
OUT = os.path.join("experiments", "my_outputs", f"lake_migori_{TAG}")
CONST = 0.1635
SHOW = [8.9797, 11.9865, 16.0]            # xTLH values drawn in the time series
SURVEY = {"RDT": {2015.55: (35.6, 15.9, 54.4), 2020.9: (17.4, 7.5, 28.1)},
          "microscopy": {2015.55: (31.1, 13.6, 50.9), 2020.9: (13.6, 4.7, 21.9)}}
CHANNEL = {"RDT": "PCR_Parasite_Prevalence", "microscopy": "Blood_Smear_Parasite_Prevalence"}
INK, MUTE, GRID = "#1f2933", "#6b7785", "#e3e7eb"
SHADES = {8.9797: "#ef8a8a", 11.9865: "#d7191c", 16.0: "#7a0a0c"}
plt.rcParams.update({"font.size": 10, "text.color": INK, "axes.labelcolor": INK, "xtick.color": MUTE, "ytick.color": MUTE})

sp = pd.read_csv(os.path.join(OUT, "SpatialReportMalariaFiltered_children.csv"))
sp = sp[sp["node"] == 3].copy()
sp["xTLH"] = sp["xTLH"].astype(float).round(4)
J = pd.read_csv(f"J_per_seed_{TAG}.csv"); J["xTLH"] = J["xTLH"].astype(float).round(4)
Y = pd.read_csv(f"scores_per_seed_year_{TAG}.csv"); Y["xTLH"] = Y["xTLH"].astype(float).round(4)

def clean(ax):
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(MUTE)

def series_panel(ax, test, title):
    clean(ax)
    ch = CHANNEL[test]
    for x in SHOW:
        t = sp[sp["xTLH"] == x]
        if t.empty: continue
        w = t.pivot_table(index="time", columns="Run_Number", values=ch) / CONST * 100
        w = w.rolling(30, center=True, min_periods=10).mean()
        yr = 2010 + (w.index.values - 1) / 365
        ax.fill_between(yr, w.min(axis=1), w.max(axis=1), color=SHADES[x], alpha=0.15, lw=0)
        ax.plot(yr, w.median(axis=1), color=SHADES[x], lw=1.8, label=f"xTLH {x:.1f}")
    for yr, (p, lo, hi) in SURVEY[test].items():
        ax.errorbar([yr], [p], yerr=[[p - lo], [hi - p]], fmt="o", color="black", ms=8, capsize=4, lw=1.8, zorder=5)
    ax.plot([], [], "o", color="black", label=f"{test} survey (95% CI)")
    ax.set_xlim(2010, 2021); ax.set_ylim(0, 70); ax.set_xlabel("Year"); ax.set_ylabel("Prevalence, children 6-59 months (%)")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")

fig, axs = plt.subplots(2, 2, figsize=(13, 9), dpi=200)
series_panel(axs[0, 0], "RDT", "A. Model PCR through time vs RDT survey")
series_panel(axs[0, 1], "microscopy", "B. Model blood-smear through time vs microscopy survey")

ax = axs[1, 0]; clean(ax)
for test, ls, mk, c in (("RDT", "-", "o", "black"), ("microscopy", "--", "s", "#d7191c")):
    t = J[J["test"] == test].groupby("xTLH")["J"].agg(["median", "min", "max"])
    ax.fill_between(t.index, t["min"], t["max"], color=c, alpha=0.12, lw=0)
    ax.plot(t.index, t["median"], color=c, ls=ls, marker=mk, ms=5, lw=2, label=f"vs {test}: median, min-max of seeds")
    b = t["median"].idxmin(); ax.plot([b], [t.loc[b, "median"]], "o", ms=12, mfc="none", mec=INK, mew=1.4)
ax.axhline(4, color="#9aa5b1", lw=1, ls=":"); ax.text(0.52, 4.4, "J = 4: close fit", color=MUTE, fontsize=8.5)
ax.set_xscale("log"); ax.set_xticks([0.5, 1, 2, 4, 8, 16]); ax.set_xticklabels(["0.5", "1", "2", "4", "8", "16"]); ax.minorticks_off()
ax.set_ylim(-1.5, 26); ax.set_xlabel("x_Temporary_Larval_Habitat (xTLH)"); ax.set_ylabel("J (lower = better fit)")
ax.set_title("C. Fit to both surveys (rings = lowest median J)", loc="left", fontsize=11, fontweight="bold")
ax.legend(frameon=False, fontsize=8.5, loc="upper right")

ax = axs[1, 1]; clean(ax)
r = Y[Y["test"] == "RDT"]
for i, x in enumerate(SHOW):
    for k, (yr, col) in enumerate(((2015, "#d7191c"), (2020, "black"))):
        d = r[(r["xTLH"] == x) & (r["year"] == yr)]
        xs = i + (k - 0.5) * 0.3 + np.linspace(-0.05, 0.05, max(len(d), 1))[: len(d)]
        ax.plot(xs, d["p_model"] * 100, "o", color=col, ms=6, alpha=0.8, label=f"model {yr}" if i == 0 else None)
for k, (yr, col) in enumerate(((2015, "#d7191c"), (2020, "black"))):
    p, lo, hi = SURVEY["RDT"][2015.55 if yr == 2015 else 2020.9]
    ax.axhspan(lo, hi, color=col, alpha=0.07); ax.axhline(p, color=col, ls="--", lw=1.4); ax.axhline(lo, color=col, ls=":", lw=0.8); ax.axhline(hi, color=col, ls=":", lw=0.8)
    ax.text(-0.58, p + 0.8, f"{yr} survey {p}%", color=col, fontsize=8.5, ha="left")
ax.set_xticks(range(len(SHOW))); ax.set_xticklabels([f"xTLH {x:.1f}" for x in SHOW])
ax.set_xlim(-0.6, 2.6); ax.set_ylim(0, 68); ax.set_ylabel("Model PCR prevalence in the survey window (%)")
ax.set_title("D. Each seed vs the RDT survey (shaded: 95% CI)", loc="left", fontsize=11, fontweight="bold")
ax.legend(frameon=False, fontsize=8.5, loc="upper left")

fig.suptitle(f"Migori County (node 3): EMOD calibration to the Kenya MIS 2015 and 2020 [{TAG}]", fontsize=13, fontweight="bold", x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.96))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"Migori_calibration_{TAG}.{ext}"), bbox_inches="tight")
print("Saved:", os.path.join(OUT, f"Migori_calibration_{TAG}.png"))
