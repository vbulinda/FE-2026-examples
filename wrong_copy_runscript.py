"""
analyzer_pickup2.py  -  analysis of pickup v2 (30-year burn-in, then 2010-2020 with ITN + case management from MIS tables)

Run from FE-2026-examples (needs analyzer_W2.py and analyzer_spatial.py there; imports only the classes):
    python analyzer_pickup2.py

Writes to <job_directory>/my_outputs/lake_pickup2/ :
  All_Age_InsetChart.csv, PfPR_ClinicalIncidence_monthly.csv, SpatialReportMalariaFiltered_all_ages.csv
  pickup_by_xTLH_county.csv     per county and xTLH: prevalence and clinical incidence, last 12 months
  pickup_prevalence_by_xTLH.png prevalence over time, one panel per county, one line per xTLH
  pickup_children_by_xTLH.csv   children under 5 vs all ages, last 12 months, whole zone
"""
import os, glob
import pandas as pd
import numpy as np
from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform

import manifest
from analyzer_W2 import InsetChartAnalyzer, MonthlyPfPRAnalyzer
from analyzer_spatial import SpatialAnalyzer

EXP_ID = "202b3b52-ea69-4b41-a270-dfb6505687b6"     # pickup v2 (ITN + case management, MIS tables)
NAME = "lake_pickup2"
START_YEAR = 2010
NODE_NAMES = {1: "Kisumu", 2: "Siaya", 3: "Migori", 4: "Homa Bay",
              5: "Kakamega", 6: "Busia", 7: "Bungoma", 8: "Vihiga"}

wdir = os.path.join(manifest.job_directory, "my_outputs", NAME)
os.makedirs(wdir, exist_ok=True)
sweep = ["Run_Number", "xTLH"]
channels = ["Statistical Population", "True Prevalence", "New Clinical Cases",
            "Infectious Vectors", "Rainfall", "Air Temperature"]
spatial_channels = ["Population", "PCR_Parasite_Prevalence", "New_Clinical_Cases"]

# does the pickup contain the spatial files?
have_spatial = False
for ed in [d for d in glob.glob(os.path.join(manifest.job_directory, "**", "*" + EXP_ID + "*"), recursive=True) if os.path.isdir(d)][:1]:
    outs = glob.glob(os.path.join(ed, "*", "output"))
    if outs:
        have_spatial = all(os.path.exists(os.path.join(outs[0], f"SpatialReportMalariaFiltered_all_ages_{c}.bin")) for c in spatial_channels)
        print("Output files of first simulation:", sorted(os.listdir(outs[0]))[:12], "...")
print("Spatial files present:", have_spatial)

with Platform("SLURM_LOCAL", job_directory=manifest.job_directory) as platform:
    analyzers = [InsetChartAnalyzer(expt_name=NAME, channels=channels, sweep_variables=sweep,
                                    start_year=START_YEAR, working_dir=wdir),
                 MonthlyPfPRAnalyzer(expt_name=NAME, sweep_variables=sweep,
                                     start_year=START_YEAR, working_dir=wdir)]
    if have_spatial:
        analyzers.append(SpatialAnalyzer(dir_name=NAME, f_base="SpatialReportMalariaFiltered", f_suffix="all_ages",
                                         exp_id=EXP_ID, spatial_channels=spatial_channels,
                                         sweep_variables=sweep, working_dir=wdir))
    AnalyzeManager(configuration={}, ids=[(EXP_ID, ItemType.EXPERIMENT)],
                   analyzers=analyzers, partial_analyze_ok=True).analyze()

if have_spatial:
    sp = pd.read_csv(os.path.join(wdir, "SpatialReportMalariaFiltered_all_ages.csv"))
    sp["xTLH"] = sp["xTLH"].astype(float).round(3)
    tmax = sp["time"].max()
    last = sp[sp["time"] > tmax - 365]                      # last 12 months of the pickup
    # per simulation first (xTLH, node, seed), then average over seeds, so seeds are not added together
    per_sim = last.groupby(["xTLH", "node", "Run_Number"]).agg(
        pcr_prevalence=("PCR_Parasite_Prevalence", "mean"),
        population=("Population", "mean"),
        cases_year=("New_Clinical_Cases", "sum")).reset_index()
    tab = per_sim.groupby(["xTLH", "node"])[["pcr_prevalence", "population", "cases_year"]].mean().reset_index()
    tab["clinical_incidence_per_1000_year"] = tab["cases_year"] / tab["population"] * 1000
    tab["county"] = tab["node"].map(NODE_NAMES)
    tab = tab[["xTLH", "node", "county", "pcr_prevalence", "clinical_incidence_per_1000_year", "population"]].round(4)
    tab.to_csv(os.path.join(wdir, "pickup_by_xTLH_county.csv"), index=False)

    print("\nMean PCR prevalence, last 12 months, by county (rows) and xTLH (columns):")
    print(tab.pivot(index="county", columns="xTLH", values="pcr_prevalence").round(3))
    print("\nClinical cases per 1,000 people per year, last 12 months:")
    print(tab.pivot(index="county", columns="xTLH", values="clinical_incidence_per_1000_year").round(0))
    print("\nAll counties together (mean of counties):")
    print(tab.groupby("xTLH")[["pcr_prevalence", "clinical_incidence_per_1000_year"]].mean().round(3))

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 4, figsize=(14, 6), sharex=True, sharey=True)
        xs = sorted(sp["xTLH"].unique())
        col = {x: plt.get_cmap("tab20")(i % 20) for i, x in enumerate(xs)}
        first_run = sp["Run_Number"].min()
        for ax, (node, gg) in zip(axs.ravel(), sp.groupby("node")):
            for (x, run), g2 in gg.groupby(["xTLH", "Run_Number"]):
                g2 = g2.sort_values("time")
                ax.plot(g2["time"] / 365, g2["PCR_Parasite_Prevalence"].rolling(14, min_periods=1).mean(),
                        lw=1, color=col[x], label=f"xTLH {x}" if run == first_run else None)
            ax.set_title(NODE_NAMES.get(node, node))
        axs[0, 0].legend(fontsize=7)
        fig.supxlabel("years since pickup start"); fig.supylabel("PCR prevalence (14-day mean)")
        fig.tight_layout(); p = os.path.join(wdir, "pickup_prevalence_by_xTLH.png"); fig.savefig(p, dpi=120)
        print("\nPlot:", p)
    except Exception as e:
        print("(plot skipped:", e, ")")
# ---------- children under 5 (closest output to the 6-59 month calibration target) ----------
mp = os.path.join(wdir, NAME, "PfPR_ClinicalIncidence_monthly.csv")
if os.path.exists(mp):
    m = pd.read_csv(mp)
    m["xTLH"] = m["xTLH"].astype(float).round(3)
    bins = sorted(m["agebin"].unique())
    print("\nAge bins in the summary report:", bins)
    u5 = bins[1] if len(bins) > 1 else bins[0]               # second bin = 0.25 to 5 years
    mmax = m["month"].max()
    L = m[m["month"] > mmax - 12]
    ps = L.groupby(["xTLH", "Run_Number", "agebin"]).agg(PfPR=("PfPR", "mean"), incidence=("Cases", "mean")).reset_index()
    ps = ps.groupby(["xTLH", "agebin"])[["PfPR", "incidence"]].mean().reset_index()
    ps["incidence_per_1000_year"] = ps["incidence"] * 1000
    kids = ps[ps["agebin"] == u5].set_index("xTLH")[["PfPR", "incidence_per_1000_year"]].add_suffix("_under5")
    allp = ps[ps["agebin"] == bins[-1]].set_index("xTLH")[["PfPR", "incidence_per_1000_year"]].add_suffix("_all_ages")
    out = kids.join(allp).round(3)
    out.to_csv(os.path.join(wdir, "pickup_children_by_xTLH.csv"))
    print("\nWhole zone, last 12 months (age bin", u5, "= children under 5, 3-59 months; bin", bins[-1], "= all ages):")
    print(out)
else:
    print("\n(no monthly summary file found: skipped the children table)")

print("\nFiles in", wdir, ":", sorted(os.listdir(wdir)))