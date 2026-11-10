"""
analyzer_pickup.py  -  analysis of the Lake Victoria pickup test (Stage 9)

Run from FE-2026-examples (needs analyzer_W2.py and analyzer_spatial.py there; imports only the classes):
    python analyzer_pickup.py

The pickup has one simulation per larval-habitat value (xTLH), so results are grouped by xTLH.
Writes to <job_directory>/my_outputs/lake_pickup30/ :
  All_Age_InsetChart.csv, PfPR_ClinicalIncidence_monthly.csv, SpatialReportMalariaFiltered_all_ages.csv
  pickup_by_xTLH_county.csv     per county and xTLH: prevalence and clinical incidence, last 12 months
  pickup_prevalence_by_xTLH.png prevalence over time, one panel per county, one line per xTLH
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

EXP_ID = "02eeb191-8ed2-4c4b-9b34-998ce91df0b0"     # pickup experiment
NAME = "lake_pickup30"
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
        for ax, (node, gg) in zip(axs.ravel(), sp.groupby("node")):
            for x, g2 in gg.groupby("xTLH"):
                ax.plot(g2["time"] / 365, g2["PCR_Parasite_Prevalence"].rolling(14, min_periods=1).mean(),
                        lw=1, label=f"xTLH {x}")
            ax.set_title(NODE_NAMES.get(node, node))
        axs[0, 0].legend(fontsize=7)
        fig.supxlabel("years since pickup start"); fig.supylabel("PCR prevalence (14-day mean)")
        fig.tight_layout(); p = os.path.join(wdir, "pickup_prevalence_by_xTLH.png"); fig.savefig(p, dpi=120)
        print("\nPlot:", p)
    except Exception as e:
        print("(plot skipped:", e, ")")
print("\nFiles in", wdir, ":", sorted(os.listdir(wdir)))