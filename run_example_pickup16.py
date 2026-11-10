"""
run_lake_pickup3.py  -  calibration pickup: pickup v2 (case management + ITN from the MIS tables) on the 20-simulation calibration burn-in,
with an extra per-county report for children 6-59 months

Starts from the 30-year burn-in (BURNIN_EXP_ID), carries each simulation's x_Temporary_Larval_Habitat and
seed, and runs 2010-2020 (PICKUP_YEARS = 11) with, in every node:
  - case management from interventions/case_management_v2.csv   (rows with phase == "pickup", one campaign per year)
  - ITN campaigns from interventions/ITN.csv + ITN_age.csv   (usage depends on age; kill/block/decay are PLACEHOLDERS)
  - IRS from interventions/IRS.json, only if that file exists (it does not until NMCP data arrive or you build it)
Tables with a single node (node_id 1) are copied to all 8 nodes; tables with all 8 node_ids are used as they are.

Put these files in <manifest.input_dir>/interventions/ :  case_management_v2.csv, ITN.csv, ITN_age.csv
Run from the FE-2026-examples folder with (emod-venv) active:   python run_lake_pickup3.py
"""
import os
import sys
import json
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
import emodpy_malaria.interventions.usage_dependent_bednet as itn
import emodpy_malaria.interventions.irs as irs
from emodpy_malaria.reporters.builtin import *

import manifest
sys.path.append(os.getcwd())
from utils_slurm import build_burnin_df

# ---------------- settings ----------------
BURNIN_EXP_IDS  = ["10b92b9f-30c4-4da0-b315-015bef183702",   # first calibration burn-in (14 runs have a saved state)
                   "4218060a-da00-43ec-9a38-18aea520105c", "68f6ffd4-4891-4831-8598-00ab12343759"]              # ID printed by run_lake_burnin_cal_rerun.py (the missing runs)
SERIALIZE_YEARS = 30         # must equal SERIALIZE_YEARS used in the burn-in
PICKUP_YEARS    = 11         # 2010-2020: covers the ITN rounds of 2014, 2015, 2020 (day 1 = 1 Jan 2010)
TEST_POP        = 5000       # must equal TEST_POP used in run_lake_burnin_cal.py
ID_REF          = "LAKE_VICTORIA"
CLIMATE_DIR     = os.path.join("climate", "LAKE", "2010001-2019365")
NODES_FILE      = os.path.join(manifest.input_dir, "demographics", "lake_nodes.csv")
INTERV_DIR      = os.path.join(manifest.input_dir, "interventions")
SPECIES         = ["gambiae", "arabiensis", "funestus"]
EVENTS_TO_RECORD = ["Received_Treatment", "Bednet_Got_New_One", "Bednet_Using", "Bednet_Discarded"]
FILES = {"air": "dtk_15arcmin_air_temperature_daily", "land": "dtk_15arcmin_land_temperature_daily",
         "rh": "dtk_15arcmin_relative_humidity_daily", "rain": "dtk_15arcmin_rainfall_daily"}
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
    config.parameters.Simulation_Duration = PICKUP_YEARS * 365 + 30   # +30 lets the last monthly report finish
    config.parameters.Custom_Individual_Events = EVENTS_TO_RECORD

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


def _num(v):
    """float() that also reads 'PLACEHOLDER_<number>' tags (unknown ITN efficacy values)."""
    if isinstance(v, str) and v.upper().startswith("PLACEHOLDER_"):
        if not getattr(_num, "warned", False):
            print("WARNING: PLACEHOLDER values are used for ITN kill/block/decay. Replace with cited values.")
            _num.warned = True
        return float(v.split("_", 1)[1])
    return float(v)


def _node_table(df, node_ids):
    """A table with a single node is copied to every node."""
    ids = set(int(n) for n in df["node_id"].unique())
    if len(ids) == 1 and ids != set(node_ids):
        print(f"NOTE: table has only node_id={sorted(ids)}; copying its rows to all {len(node_ids)} nodes")
        df = pd.concat([df.assign(node_id=int(n)) for n in node_ids], ignore_index=True)
    return df


def build_camp():
    camp.set_schema(manifest.schema_file)
    node_ids = pd.read_csv(NODES_FILE)["node_id"].tolist()

    # ---- case management (one treatment-seeking campaign per node and year) ----
    cm_all = _node_table(pd.read_csv(os.path.join(INTERV_DIR, "case_management_v2.csv")), node_ids)
    for node in cm_all["node_id"].unique():
        n_df = cm_all[(cm_all["node_id"] == node) & (cm_all["phase"] == "pickup")].reset_index(drop=True)
        for year in n_df["year"].unique():
            sub = n_df[n_df["year"] == year].reset_index(drop=True)
            targets = [{"trigger": str(sub["trigger"][r]), "coverage": float(sub["coverage"][r]),
                        "agemin": float(sub["age_min"][r]), "agemax": float(sub["age_max"][r]),
                        "rate": float(sub["rate"][r])} for r in range(len(sub))]
            cm.add_treatment_seeking(camp, node_ids=[int(node)], start_day=int(sub["start_day"][0]),
                                     duration=int(sub["duration"][0]), drug=["Artemether", "Lumefantrine"],
                                     targets=targets, broadcast_event_name="Received_Treatment")

    # ---- ITN: usage depends on age; kill/block/decay may be placeholders ----
    itn_all = _node_table(pd.read_csv(os.path.join(INTERV_DIR, "ITN_v8.csv")), node_ids)
    itn_age = pd.read_csv(os.path.join(INTERV_DIR, "ITN_age.csv"))
    season_path = os.path.join(INTERV_DIR, "ITN_season.csv")
    if os.path.exists(season_path):
        sdf = pd.read_csv(season_path)
        seasonal = {"Times": list(sdf["season_time"]), "Values": list(sdf["season_usage"])}
    else:
        print("NOTE: ITN_season.csv not found; using flat seasonal usage.")
        seasonal = {"Times": [0, 365], "Values": [1.0, 1.0]}
    for node in itn_all["node_id"].unique():
        n_df = itn_all[(itn_all["node_id"] == node) & (itn_all["phase"] == "pickup")].reset_index(drop=True)
        for year in n_df["year"].unique():
            sub = n_df[n_df["year"] == year].reset_index(drop=True)
            age_y = itn_age[itn_age["year"] == year]
            if "node_id" in itn_age.columns:      # county-specific usage by age
                age_y = age_y[age_y["node_id"] == int(node)]
            discard = {"Expiration_Period_Distribution": "WEIBULL_DISTRIBUTION",
                       "Expiration_Period_Kappa": _num(sub["discard_k"][0]),
                       "Expiration_Period_Lambda": _num(sub["discard_l"][0])}
            itn.add_scheduled_usage_dependent_bednet(
                camp, node_ids=[int(node)], intervention_name="UsageDependentBednet",
                start_day=int(sub["start_day"][0]), demographic_coverage=float(sub["coverage"][0]),
                killing_initial_effect=_num(sub["kill_effect"][0]),
                killing_decay_time_constant=int(_num(sub["kill_decay"][0])),
                blocking_initial_effect=_num(sub["block_effect"][0]),
                blocking_decay_time_constant=int(_num(sub["block_decay"][0])),
                age_dependence={"Times": list(age_y["age"]), "Values": list(age_y["age_usage"])},
                seasonal_dependence=seasonal, discard_config=discard)

    # ---- IRS: only if an IRS.json exists (NMCP data or a tagged placeholder you built) ----
    irs_path = os.path.join(INTERV_DIR, "IRS.json")
    if os.path.exists(irs_path):
        with open(irs_path) as f:
            irs_all = _node_table(pd.DataFrame(json.load(f)), node_ids)
        for node in irs_all["node_id"].unique():
            n_df = irs_all[(irs_all["node_id"] == node) & (irs_all["phase"] == "pickup")].reset_index(drop=True)
            for r in range(len(n_df)):
                irs.add_scheduled_irs_housing_modification(
                    camp, node_ids=[int(node)], start_day=int(n_df["start_day"][r]),
                    demographic_coverage=float(n_df["coverage"][r]),
                    killing_initial_effect=float(n_df["kill_effect"][r]), killing_box_duration=int(n_df["kill_duration"][r]),
                    killing_decay_time_constant=int(n_df["kill_decay"][r]),
                    repelling_initial_effect=float(n_df["repel_effect"][r]), repelling_box_duration=int(n_df["repel_duration"][r]),
                    repelling_decay_time_constant=int(n_df["repel_decay"][r]))
        print("IRS included from", irs_path)
    else:
        print("NOTE: no IRS.json found: this run has NO IRS.")

    fix_property_restrictions(camp.campaign_dict)      # null Node_Property_Restrictions -> []
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
                        partition=manifest.partition, time="12:00:00",
                        modules=[manifest.singularity_module],
                        max_running_jobs=manifest.max_running_jobs)

    task = EMODTask.from_default2(
        config_path="config.json", eradication_path=manifest.eradication_path,
        campaign_builder=build_camp, schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn, ep4_custom_cb=None,
        demog_builder=build_demog, plugin_report=None)
    task.set_sif(manifest.SIF_PATH, platform)
    task.common_assets.add_directory(os.path.join(manifest.input_dir, "climate"), relative_path="climate")

    need = [os.path.join(INTERV_DIR, n) for n in ("case_management_v2.csv", "ITN_v8.csv", "ITN_age.csv")]
    miss = [p for p in need if not os.path.exists(p)]
    if miss:
        raise SystemExit("Missing intervention files:\n  " + "\n  ".join(miss) + "\nUpload them to " + INTERV_DIR)

    # one pickup simulation per burn-in simulation
    if any("PASTE" in i for i in BURNIN_EXP_IDS):
        raise SystemExit("Set BURNIN_EXP_IDS: paste the ID printed by run_lake_burnin_cal_rerun.py, then run again.")
    parts = [build_burnin_df(i, platform, SERIALIZE_YEARS * 365) for i in BURNIN_EXP_IDS]
    burnin_df = pd.concat(parts, ignore_index=True)
    # keep only simulations that actually have their saved state, one per (habitat value, seed)
    has = [all(os.path.exists(os.path.join(burnin_df["serialized_file_path"][i], "output", fn))
               for fn in burnin_df["Serialized_Population_Filenames"][i]) for i in burnin_df.index]
    burnin_df = burnin_df[has].copy()
    burnin_df["x_Temporary_Larval_Habitat"] = burnin_df["x_Temporary_Larval_Habitat"].astype(float).round(4)
    burnin_df = burnin_df.drop_duplicates(subset=["x_Temporary_Larval_Habitat", "Run_Number"]).reset_index(drop=True)
    if len(burnin_df) != 20:
        print("WARNING: expected 20 burn-in states, found", len(burnin_df), "- check the burn-in IDs before going on.")
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
                                 event_trigger_list=EVENTS_TO_RECORD,
                                 filename_suffix="_".join(("node", str(n))))
    add_malaria_summary_report(task, manifest, start_day=1, end_day=end_day, reporting_interval=30,
                               age_bins=[0.25, 5, 115], max_number_reports=PICKUP_YEARS * 13 + 1,
                               filename_suffix="monthly", pretty_format=True)
    add_spatial_report_malaria_filtered(task, manifest, start_day=1, end_day=end_day, reporting_interval=1,
                                        node_ids=None, min_age_years=0.25, max_age_years=100,
                                        spatial_output_channels=["Population", "PCR_Parasite_Prevalence",
                                                                 "New_Clinical_Cases"],
                                        filename_suffix="all_ages")
    # children 6-59 months, per county: the calibration target (MIS prevalence is measured in this age group)
    add_spatial_report_malaria_filtered(task, manifest, start_day=1, end_day=end_day, reporting_interval=1,
                                        node_ids=None, min_age_years=0.5, max_age_years=5,
                                        spatial_output_channels=["Population", "PCR_Parasite_Prevalence",
                                                                 "New_Clinical_Cases"],
                                        filename_suffix="children")

    user = os.getlogin()
    experiment = Experiment.from_builder(builder, task, name=f"{user}_LAKE_pickup_cal14")
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