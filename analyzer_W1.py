"""
Stage 4 analysis: turns the outputs of run_outputs.py into CSV files.
Needs analyzer_W2.py (from solution_scripts/Week2) and analyzer_spatial.py (repo root)
in the same folder. Set EXP_ID to the experiment id printed by run_outputs.py.
"""
import os
import pandas as pd
from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform

import manifest
from analyzer_W2 import InsetChartAnalyzer, MonthlyPfPRAnalyzer
from analyzer_spatial import SpatialAnalyzer

EXP_ID = "b057ec0e-d1b8-4299-8ecf-e3e23f2015f9"
NAME = "lake_outputs_test"
START_YEAR = 2010

wdir = os.path.join(manifest.job_directory, "my_outputs", NAME)
os.makedirs(wdir, exist_ok=True)
sweep = ["Run_Number"]

# ---- look at what the simulation actually wrote -----------------------------
import glob
spatial_channels = ["Population", "PCR_Parasite_Prevalence", "New_Clinical_Cases"]
have_spatial = False
exp_dirs = [d for d in glob.glob(os.path.join(manifest.job_directory, "**", EXP_ID), recursive=True) if os.path.isdir(d)]
print("Experiment folder(s):", exp_dirs)
for ed in exp_dirs[:1]:
    outs = glob.glob(os.path.join(ed, "*", "output"))
    if outs:
        files = sorted(os.listdir(outs[0]))
        print("Output files of first simulation:\n ", "\n  ".join(files))
        have_spatial = all(os.path.exists(os.path.join(outs[0], f"SpatialReportMalariaFiltered_all_ages_{c}.bin")) for c in spatial_channels)
print("Spatial files present:", have_spatial)
channels = ["Statistical Population", "True Prevalence", "New Clinical Cases",
            "Infectious Vectors", "Rainfall", "Air Temperature"]

with Platform("SLURM_LOCAL", job_directory=manifest.job_directory) as platform:
    analyzers = [
        InsetChartAnalyzer(expt_name=NAME, channels=channels, sweep_variables=sweep,
                           start_year=START_YEAR, working_dir=wdir),
        MonthlyPfPRAnalyzer(expt_name=NAME, sweep_variables=sweep,
                            start_year=START_YEAR, working_dir=wdir),
    ]
    if have_spatial:
        analyzers.append(SpatialAnalyzer(dir_name=NAME, f_base="SpatialReportMalariaFiltered", f_suffix="all_ages",
                        exp_id=EXP_ID, spatial_channels=spatial_channels,
                        sweep_variables=sweep, working_dir=wdir))
    else:
        print("\nNo spatial files in this experiment: skipping SpatialAnalyzer. Re-run run_outputs.py to get them.")
    AnalyzeManager(configuration={}, ids=[(EXP_ID, ItemType.EXPERIMENT)],
                   analyzers=analyzers, partial_analyze_ok=True).analyze()

if have_spatial:
    sp = pd.read_csv(os.path.join(wdir, "SpatialReportMalariaFiltered_all_ages.csv"))
    print("\nSpatial rows:", len(sp), " nodes:", sorted(sp["node"].unique()))
    print(sp.groupby("node")[["Population", "PCR_Parasite_Prevalence", "New_Clinical_Cases"]].mean().round(3))
print("\nFiles in", wdir, ":", os.listdir(wdir))