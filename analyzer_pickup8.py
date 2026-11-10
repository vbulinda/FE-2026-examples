"""
analyzer_pickup3.py - collect the calibration pickup (run_example_pickup3.py) into CSV files.
Run from FE-2026-examples (needs analyzer_W2.py and analyzer_spatial.py there):
    python analyzer_pickup3.py
Writes to <job_directory>/my_outputs/lake_pickup_cal/ :
    SpatialReportMalariaFiltered_all_ages.csv   per node, day, seed, xTLH (all ages)
    SpatialReportMalariaFiltered_children.csv   same, children 6-59 months (the calibration target)
"""
import os, glob
import pandas as pd
from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform
import manifest
from analyzer_spatial import SpatialAnalyzer

EXP_ID = "7ff5eca3-7b03-476c-886c-8458df5b241e"
NAME = "lake_pickup_cal6"
wdir = os.path.join(manifest.job_directory, "my_outputs", NAME)
os.makedirs(wdir, exist_ok=True)
sweep = ["Run_Number", "xTLH"]
spatial_channels = ["Population", "PCR_Parasite_Prevalence", "New_Clinical_Cases"]

# which spatial files exist in the first finished simulation?
suffixes = []
for ed in [d for d in glob.glob(os.path.join(manifest.job_directory, "**", "*" + EXP_ID + "*"), recursive=True) if os.path.isdir(d)][:1]:
    outs = sorted(glob.glob(os.path.join(ed, "*", "output")))
    for o in outs:
        names = os.listdir(o)
        if any(n.startswith("SpatialReportMalariaFiltered_children") for n in names):
            print("Output files of one simulation:", sorted(names)[:14]); break
    for s in ("all_ages", "children"):
        if outs and all(os.path.exists(os.path.join(outs[0], f"SpatialReportMalariaFiltered_{s}_{c}.bin")) for c in spatial_channels):
            suffixes.append(s)
print("Spatial report sets found:", suffixes)
if "children" not in suffixes:
    raise SystemExit("The children report is missing in the first simulation: the pickup may not have finished. Check status first.")

# one AnalyzeManager pass per report set: two spatial analyzers in one pass clash (KeyError on the second file set)
done = []
with Platform("SLURM_LOCAL", job_directory=manifest.job_directory) as platform:
    for s in ["children"] + [x for x in suffixes if x != "children"]:
        print(f"\n=== analysing the '{s}' report ===")
        try:
            an = SpatialAnalyzer(dir_name=NAME, f_base="SpatialReportMalariaFiltered", f_suffix=s,
                                 exp_id=EXP_ID, spatial_channels=spatial_channels,
                                 sweep_variables=sweep, working_dir=wdir)
            want = [f"output/SpatialReportMalariaFiltered_{s}_{c}.bin" for c in spatial_channels]
            print("analyzer will load:", getattr(an, "filenames", "(no .filenames attribute)"))
            if sorted(getattr(an, "filenames", [])) != sorted(want):
                print("-> correcting the file list to:", want)
                an.filenames = want
            AnalyzeManager(configuration={}, ids=[(EXP_ID, ItemType.EXPERIMENT)],
                           analyzers=[an], partial_analyze_ok=True).analyze()
            done.append(s)
        except Exception as e:
            print(f"!! '{s}' failed: {type(e).__name__}: {e}")
print("Report sets analysed:", done)
if "children" not in done:
    raise SystemExit("The children analysis failed. Send me the lines above.")

print("\nFiles written:")
for f in sorted(glob.glob(os.path.join(wdir, "**", "SpatialReportMalariaFiltered_*.csv"), recursive=True)):
    d = pd.read_csv(f, nrows=3); print(" ", f, list(d.columns))