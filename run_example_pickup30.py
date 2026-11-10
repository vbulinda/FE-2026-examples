"""
run_lake_pickup30.py  -  Stage 8: pickup from the REAL 30-year burn-in, Lake Victoria, 8 counties

Starts from the saved population states of a burn-in (BURNIN_EXP_ID), carries over that burn-in's
x_Temporary_Larval_Habitat and seed, then runs PICKUP_YEARS more years with:
  - case management (treatment seeking) in all 8 nodes   <- PLACEHOLDER coverage, see CM_* below
  - the same reports as before (monthly summary, births per node, daily spatial report)
ITN and IRS are NOT added yet. Add them in version 2 once this version succeeds.

One pickup simulation is made per burn-in simulation (24 for the real burn-in).
Run from the FE-2026-examples folder (needs manifest.py and utils_slurm.py there):
    python run_lake_pickup.py
"""
import os
import sys
import pathlib
from functools import partial
import pandas as pd

from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from emodpy.emod_task import EMODTask
import emod_api.campaign as camp
import emodpy_malaria.demographics.MalariaDemographics as Demographics
import emod_api.demographics.PreDefinedDistributions as Distributions
import emodpy_malaria.interventions.treatment_seeking as cm
from emodpy_malaria.reporters.builtin import *

import manifest
sys.path.append(os.getcwd())
from utils_slurm import build_burnin_df

# ---------------- settings ----------------
BURNIN_EXP_ID   = "2d3e47b4-199b-46d5-a19f-6cc14bd929a9"   # real 30-year burn-in (24 simulations)
SERIALIZE_YEARS = 30         # must equal SERIALIZE_YEARS used in the burn-in
PICKUP_YEARS    = 5          # years simulated after the burn-in (2010 climate onward)
TEST_POP        = 1000       # must equal TEST_POP used in the burn-in
ID_REF          = "LAKE_VICTORIA"
CLIMATE_DIR     = os.path.join("climate", "LAKE", "2010001-2019365")
NODES_FILE      = os.path.join(manifest.input_dir, "demographics", "lake_nodes.csv")
SPECIES         = ["gambiae", "arabiensis", "funestus"]
FILES = {"air": "dtk_15arcmin_air_temperature_daily", "land": "dtk_15arcmin_land_temperature_daily",
         "rh": "dtk_15arcmin_relative_humidity_daily", "rain": "dtk_15arcmin_rainfall_daily"}
# PLACEHOLDER care-seeking coverage (fraction of cases treated). Replace with Lake-zone values from the MIS.
CM_U5, CM_5_15, CM_15PLUS, CM_SEVERE = 0.40, 0.24, 0.16, 0.80
# ------------------------------------------


def set_param_fn(config):
    import emodpy_malaria.malaria_config as conf
    config = conf.set_team_defaults(config, manifest)
    config.parameters.Climate_Model = "CLIMATE_BY_DATA"
    config.parameters.Enable_Default_Reporting = 1

    config.parameters.Air_Temperature_Filename   = os.path.join(CLIMATE_DIR, FILES["air"] + ".bin")
    config.parameters.Land_Temperature_Filename  = os.path.join(CLIMATE_DIR, FILES["land"] + ".bin")
    config.parameters.Rainfall_Filename          = os.path.join(CLIMATE_DIR, FILES["rain"] + ".bin")
    config.parameters.Relative_Humidity_Filename = os.path.join(CLIMATE_DIR, FILES["rh"] + ".bin")

    conf.add_species(config, manifest, SPECIES)
    config.parameters.Simulation_Duration = PICKUP_YEARS * 365

    # READ the burn-in state (the specific file and path are set per simulation below)
    config.parameters.Serialized_Population_Reading_Type = "READ"
    config.parameters.Serialization_Mask_Node_Read = 0
    config.parameters.Serialization_Time_Steps = [SERIALIZE_YEARS * 365]
    return config


def fix_property_restrictions(obj):
    """EMOD needs these fields to be arrays. This emodpy-malaria build writes null in some events, so
    replace anything that is not a list with an empty list (same fix as in site_setup_Kenya_LakeVictoria.py)."""
    keys = ("Node_Property_Restrictions", "Property_Restrictions", "Property_Restrictions_Within_Node")
    if isinstance(obj, dict):
        for key, value in list(obj.items()):
            if key in keys:
                if not isinstance(value, list):
                    obj[key] = []
            else:
                fix_property_restrictions(value)
    elif isinstance(obj, list):
        for item in obj:
            fix_property_restrictions(item)


def build_camp():
    camp.set_schema(manifest.schema_file)
    node_ids = pd.read_csv(NODES_FILE)["node_id"].tolist()
    targets = [{"trigger": "NewClinicalCase", "coverage": CM_U5,      "agemin": 0,  "agemax": 5,   "rate": 0.3},
               {"trigger": "NewClinicalCase", "coverage": CM_5_15,    "agemin": 5,  "agemax": 15,  "rate": 0.3},
               {"trigger": "NewClinicalCase", "coverage": CM_15PLUS,  "agemin": 15, "agemax": 115, "rate": 0.3},
               {"trigger": "NewSevereCase",   "coverage": CM_SEVERE,   "agemin": 0,  "agemax": 115, "rate": 0.5}]
    cm.add_treatment_seeking(camp, start_day=1, drug=["Artemether", "Lumefantrine"],
                             targets=targets, node_ids=node_ids,
                             broadcast_event_name="Received_Treatment")
    fix_property_restrictions(camp.campaign_dict)      # fixes the null Node_Property_Restrictions error
    return camp


def update_serialize_parameters(simulation, df, x: int):
    """Point this pickup simulation at the matching burn-in state; carry over its larval habitat and seed."""
    path = df["serialized_file_path"][x]
    xTLH = float(df["x_Temporary_Larval_Habitat"][x])
    seed = int(df["Run_Number"][x])
    simulation.task.set_parameter("Serialized_Population_Filenames", df["Serialized_Population_Filenames"][x])
    simulation.task.set_parameter("Serialized_Population_Path", os.path.join(path, "output"))
    simulation.task.set_parameter("x_Temporary_Larval_Habitat", xTLH)
    simulation.task.set_parameter("Run_Number", seed)
    return {"xTLH": xTLH, "Run_Number": seed}


def build_demog():
    """Same 8 nodes and population as the burn-in (must match for the state to load)."""
    nodes = pd.read_csv(NODES_FILE)
    nodes["pop"] = TEST_POP
    nodes.to_csv("lake_nodes_test.csv", index=False)
    demog = Demographics.from_csv(input_file="lake_nodes_test.csv", id_ref=ID_REF,
                                  init_prev=0.01, include_biting_heterogeneity=True)
    demog.SetEquilibriumVitalDynamics()
    demog.SetAgeDistribution(Distributions.AgeDistribution_SSAfrica)
    return demog


def show_failure(exp_id):
    """Print the end of the log of the first simulation folder that has an error, so the cause is visible."""
    import glob
    dirs = [d for d in glob.glob(os.path.join(manifest.job_directory, "**", "*" + exp_id + "*"), recursive=True) if os.path.isdir(d)]
    if not dirs:
        print("Experiment folder not found under", manifest.job_directory); return
    names = ["StdErr.txt", "stderr.txt", "StdOut.txt", "stdout.txt", "stderr", "stdout", "job_status.txt"]
    shown = 0
    for sim in sorted(glob.glob(os.path.join(dirs[0], "*", ""))):
        for n in names:
            p = os.path.join(sim, n)
            if os.path.exists(p) and os.path.getsize(p) > 0:
                lines = open(p, errors="replace").read().splitlines()
                print(f"\n--- {p} (last 25 lines) ---")
                print("\n".join(lines[-25:]))
                shown += 1
        if shown:
            break
    if not shown:
        print("No log files with content found. Folder listing:", os.listdir(dirs[0]))
    print("\nPaste the log lines above to Claude.")


def general_sim(selected_platform):
    platform = Platform(selected_platform, job_directory=manifest.job_directory,
                        partition=manifest.partition, time=manifest.sim_time,
                        modules=[manifest.singularity_module],
                        max_running_jobs=manifest.max_running_jobs)

    task = EMODTask.from_default2(
        config_path="config.json", eradication_path=manifest.eradication_path,
        campaign_builder=build_camp, schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn, ep4_custom_cb=None,
        demog_builder=build_demog, plugin_report=None)
    task.set_sif(manifest.SIF_PATH, platform)
    task.common_assets.add_directory(os.path.join(manifest.input_dir, "climate"), relative_path="climate")

    # one pickup simulation per burn-in simulation
    burnin_df = build_burnin_df(BURNIN_EXP_ID, platform, SERIALIZE_YEARS * 365)
    print("Burn-in simulations found:", len(burnin_df))
    print(burnin_df[["simid", "Run_Number", "x_Temporary_Larval_Habitat", "Serialized_Population_Filenames"]])
    # PREFLIGHT: the saved state must exist where the pickup will look for it
    bad = 0
    for i in burnin_df.index:
        outdir = os.path.join(burnin_df["serialized_file_path"][i], "output")
        for fn in burnin_df["Serialized_Population_Filenames"][i]:
            if not os.path.exists(os.path.join(outdir, fn)):
                bad += 1
                have = sorted(f for f in os.listdir(outdir)) if os.path.isdir(outdir) else ["(no output folder)"]
                print(f"MISSING state file: {os.path.join(outdir, fn)}\n  files there: {have}")
    if bad:
        raise SystemExit("Saved states not found, so the pickup cannot start. Paste the lines above to Claude. "
                         "Check that SERIALIZE_YEARS matches the burn-in (the file name holds the day number).")
    print("Preflight OK: all burn-in state files found.")

    builder = SimulationBuilder()
    builder.add_sweep_definition(partial(update_serialize_parameters, df=burnin_df), range(len(burnin_df.index)))

    end_day = PICKUP_YEARS * 365
    for n in pd.read_csv(NODES_FILE)["node_id"].tolist():
        add_report_event_counter(task, manifest, start_day=1, end_day=end_day, node_ids=[n],
                                 min_age_years=0, max_age_years=100,
                                 event_trigger_list=["Births"],
                                 filename_suffix="_".join(("node", str(n))))
    add_malaria_summary_report(task, manifest, start_day=1, end_day=end_day, reporting_interval=30,
                               age_bins=[0.25, 5, 115], max_number_reports=PICKUP_YEARS * 13,
                               filename_suffix="monthly", pretty_format=True)
    add_spatial_report_malaria_filtered(task, manifest, start_day=1, end_day=end_day, reporting_interval=1,
                                        node_ids=None, min_age_years=0.25, max_age_years=100,
                                        spatial_output_channels=["Population", "PCR_Parasite_Prevalence",
                                                                 "New_Clinical_Cases"],
                                        filename_suffix="all_ages")

    user = os.getlogin()
    experiment = Experiment.from_builder(builder, task, name=f"{user}_LAKE_pickup30")
    experiment.run(wait_until_done=True, platform=platform)

    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.\n")
        show_failure(str(experiment.uid))
        exit()
    print(f"Experiment {experiment.uid} succeeded with {len(burnin_df)} simulations.")
    print("Pickup ID:", experiment.uid)


if __name__ == "__main__":
    import emod_malaria.bootstrap as dtk
    dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    os.chmod(manifest.eradication_path, 0o755)
    general_sim("SLURM_LOCAL")