"""
run_bungoma_pickup_v2.py - Bungoma (node 7) pickup, 6,144 days (2010 to
31 Oct 2026), real NASA POWER weather, no IRS. Built from
run_migori_structure3.py; the changes are in the calibration log,
steps 45 and 51.   Run: python run_bungoma_pickup_v2.py
"""
import os, sys, json, pathlib
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
NODE = 7
BURNIN_EXP_IDS = ["b4f301d5-43e7-4ed8-98be-d3b6fac618a2"]
SERIALIZE_YEARS = 50
SIM_DAYS = 6144   # 1 Jan 2010 to 31 Oct 2026; equals the climate days
TEST_POP = 5000
G_MULT, ARA_MULT, FUN_MULT = 20, 0.0225, 0.006   # as in the burn-in
IRS_2017_COVERAGE = float(os.environ.get("IRS_COV", "0.81"))   # coverage used for every IRS round when IRS_REAL=1
EXPECTED_STATES = 52
NET_TAG = os.environ.get("NET_TAG", "base")
SCENARIO = os.environ.get("SCENARIO", "full")   # full, no_irs, no_itn17, no_irs_itn17
IRS_REAL = os.environ.get("IRS_REAL", "") == "1"   # 1 = IRS rounds 2017, 2018, 2019 all at IRS_2017_COVERAGE
IRS_DAYS = int(os.environ.get("IRS_DAYS", "0"))   # 0 = table value (300 days of full killing effect); else override, e.g. 150
IMPORT_PER_MONTH = int(os.environ.get("IMPORT", "1"))   # infections imported into the node every 30 days (0 = closed node)
SUBSET = os.environ.get("SUBSET", "") == "1"   # 1 = xTLH 8.9797 and 16.0, seeds 0 and 1 only
ID_REF = "LAKE_VICTORIA"
CLIMATE_DIR = os.path.join("climate", "LAKE_BUNGOMA_NASA", "2010001-2026304")
NODES_FILE = os.path.join(manifest.input_dir, "demographics", "lake_nodes_bungoma.csv")
INTERV_DIR = os.path.join(manifest.input_dir, "interventions")
SPECIES = ["gambiae", "arabiensis", "funestus"]
EVENTS_TO_RECORD = ["Received_Treatment", "Bednet_Got_New_One", "Bednet_Using", "Bednet_Discarded"]
FILES = {"air": "dtk_15arcmin_air_temperature_daily", "land": "dtk_15arcmin_land_temperature_daily",
         "rh": "dtk_15arcmin_relative_humidity_daily", "rain": "dtk_15arcmin_rainfall_daily"}
CHANNELS = ["Population", "PCR_Parasite_Prevalence", "Blood_Smear_Parasite_Prevalence", "New_Clinical_Cases"]
# ------------------------------------------


def set_param_fn(config):
    import emodpy_malaria.malaria_config as conf
    config = conf.set_team_defaults(config, manifest)
    config.parameters.Climate_Model = "CLIMATE_BY_DATA"
    config.parameters.Enable_Default_Reporting = 1
    config.parameters.Air_Temperature_Filename = os.path.join(CLIMATE_DIR, FILES["air"] + ".bin")
    config.parameters.Land_Temperature_Filename = os.path.join(CLIMATE_DIR, FILES["land"] + ".bin")
    config.parameters.Rainfall_Filename = os.path.join(CLIMATE_DIR, FILES["rain"] + ".bin")
    config.parameters.Relative_Humidity_Filename = os.path.join(CLIMATE_DIR, FILES["rh"] + ".bin")
    conf.add_species(config, manifest, SPECIES)
    sp = conf.set_max_larval_capacity
    sp(config, "gambiae", "WATER_VEGETATION", int(20000000 * G_MULT))
    sp(config, "arabiensis", "TEMPORARY_RAINFALL", int(800000000 * ARA_MULT))
    sp(config, "arabiensis", "CONSTANT", int(80000000 * ARA_MULT))
    sp(config, "funestus", "TEMPORARY_RAINFALL", int(800000000 * FUN_MULT))
    sp(config, "funestus", "CONSTANT", int(80000000 * FUN_MULT))
    config.parameters.Simulation_Duration = SIM_DAYS
    config.parameters.Custom_Individual_Events = EVENTS_TO_RECORD
    config.parameters.Serialized_Population_Reading_Type = "READ"
    config.parameters.Serialization_Mask_Node_Read = 0
    config.parameters.Serialization_Time_Steps = [SERIALIZE_YEARS * 365]
    return config


def fix_property_restrictions(obj):
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
    if isinstance(v, str) and v.upper().startswith("PLACEHOLDER_"):
        if not getattr(_num, "warned", False):
            print("WARNING: PLACEHOLDER values are used for ITN discard. Replace with cited values.")
            _num.warned = True
        return float(v.split("_", 1)[1])
    return float(v)


def _node3(df, name):
    """Keep only the rows of NODE. A table with no node 3 rows is an error, not a silent empty run."""
    out = df[df["node_id"].astype(int) == NODE].reset_index(drop=True)
    if out.empty:
        raise SystemExit(f"{name} has no rows for node_id {NODE}")
    print(f"{name}: {len(out)} rows for node {NODE}")
    return out


def build_camp():
    camp.set_schema(manifest.schema_file)
    # ---- case management ----
    cm_all = _node3(pd.read_csv(os.path.join(INTERV_DIR, "case_management_v2.csv")), "case_management_v2.csv")
    n_df = cm_all[cm_all["phase"] == "pickup"].reset_index(drop=True)
    ext = []
    for y in range(2021, 2027):   # 2020 rows carried forward, log step 51
        rows = n_df[n_df["year"] == 2020].copy()
        rows["year"] = y
        rows["start_day"] = (y - 2010) * 365 + 1
        ext.append(rows)
    n_df = pd.concat([n_df] + ext, ignore_index=True)
    for year in n_df["year"].unique():
        sub = n_df[n_df["year"] == year].reset_index(drop=True)
        targets = [{"trigger": str(sub["trigger"][r]), "coverage": float(sub["coverage"][r]),
                    "agemin": float(sub["age_min"][r]), "agemax": float(sub["age_max"][r]),
                    "rate": float(sub["rate"][r])} for r in range(len(sub))]
        cm.add_treatment_seeking(camp, node_ids=[NODE], start_day=int(sub["start_day"][0]),
                                 duration=int(sub["duration"][0]), drug=["Artemether", "Lumefantrine"],
                                 targets=targets, broadcast_event_name="Received_Treatment")
    # ---- ITN ----
    itn_all = _node3(pd.read_csv(os.path.join(INTERV_DIR, "ITN_v8.csv")), "ITN_v8.csv")
    itn_age = pd.read_csv(os.path.join(INTERV_DIR, "ITN_age.csv"))
    season_path = os.path.join(INTERV_DIR, "ITN_season.csv")
    if os.path.exists(season_path):
        sdf = pd.read_csv(season_path)
        seasonal = {"Times": list(sdf["season_time"]), "Values": list(sdf["season_usage"])}
    else:
        print("NOTE: ITN_season.csv not found; using flat seasonal usage.")
        seasonal = {"Times": [0, 365], "Values": [1.0, 1.0]}
    n_df = itn_all[itn_all["phase"] == "pickup"].reset_index(drop=True)
    base = n_df[n_df["year"] == 2017].iloc[[0]]   # PBO benefit not modelled
    NET_EXTRA = [int(v) for v in os.environ.get("NET_EXTRA", "").split(",") if v]
    n_df.loc[n_df["year"] == 2017, "start_day"] = int(os.environ.get("NET_2017", "2800"))
    extra = [(2010 + v // 365, v, 0.592) for v in NET_EXTRA]
    for y, sd, cv in [(2021, 4180, 0.825), (2024, 5291, 0.825)] + extra:
        n_df = pd.concat([n_df, base.assign(year=y, start_day=sd, coverage=cv)],
                         ignore_index=True)
    print("ITN rounds (year, start_day, coverage):")
    print(n_df[["year", "start_day", "coverage"]].to_string())
    for year in n_df["year"].unique():
        sub = n_df[n_df["year"] == year].reset_index(drop=True)
        if year == 2017 and SCENARIO in ("no_itn17", "no_irs_itn17"):
            continue
        yrs = [y for y in itn_age["year"].unique() if y <= year] or [int(itn_age["year"].min())]
        age_y = itn_age[itn_age["year"] == max(yrs)]
        if "node_id" in itn_age.columns:
            age_y = age_y[age_y["node_id"] == NODE]
        discard = {"Expiration_Period_Distribution": "WEIBULL_DISTRIBUTION",
                   "Expiration_Period_Kappa": _num(sub["discard_k"][0]),
                   "Expiration_Period_Lambda": _num(sub["discard_l"][0])}
        itn.add_scheduled_usage_dependent_bednet(
            camp, node_ids=[NODE], intervention_name="UsageDependentBednet",
            start_day=int(sub["start_day"][0]), demographic_coverage=float(sub["coverage"][0]),
            killing_initial_effect=_num(sub["kill_effect"][0]),
            killing_decay_time_constant=int(_num(sub["kill_decay"][0])),
            blocking_initial_effect=_num(sub["block_effect"][0]),
            blocking_decay_time_constant=int(_num(sub["block_decay"][0])),
            age_dependence={"Times": list(age_y["age"]), "Values": list(age_y["age_usage"])},
            seasonal_dependence=seasonal, discard_config=discard)
    # ---- IRS (2017 coverage set to IRS_2017_COVERAGE) ----
    irs_path = os.path.join(INTERV_DIR, "IRS.json")
    if False:   # Bungoma has no IRS (no node 7 rows in IRS.json)
        with open(irs_path) as f:
            irs_all = _node3(pd.DataFrame(json.load(f)), "IRS.json")
        n_df = irs_all[irs_all["phase"] == "pickup"].reset_index(drop=True)
        if IRS_REAL:   # add the 2019 (28 Jan 2019, day 3312) and 2020 (10 Feb 2020, day 3690) rounds as copies of the last row
            extra = n_df.iloc[[len(n_df) - 1]]
            n_df = pd.concat([n_df, extra.assign(start_day=3312), extra.assign(start_day=3690)], ignore_index=True)
        for r in range(len(n_df)):
            sd = int(n_df["start_day"][r]); cov = float(n_df["coverage"][r])
            if 2555 <= sd < 2920 or IRS_REAL:
                cov = IRS_2017_COVERAGE
            print(f"IRS node {NODE}: start_day {sd}, coverage used {cov} (table value {float(n_df['coverage'][r])}), full-effect days {int(IRS_DAYS or n_df['kill_duration'][r])}")
            irs.add_scheduled_irs_housing_modification(
                camp, node_ids=[NODE], start_day=sd, demographic_coverage=cov,
                killing_initial_effect=float(n_df["kill_effect"][r]), killing_box_duration=int(IRS_DAYS or n_df["kill_duration"][r]),
                killing_decay_time_constant=int(n_df["kill_decay"][r]),
                repelling_initial_effect=float(n_df["repel_effect"][r]), repelling_box_duration=int(n_df["repel_duration"][r]),
                repelling_decay_time_constant=int(n_df["repel_decay"][r]))
    else:
        print("NOTE: Bungoma has no IRS rows (node 7): this run has NO IRS.")
    if IMPORT_PER_MONTH > 0:   # imported infections: N people infected every 30 days, immunity respected
        from emodpy_malaria.interventions.outbreak import add_outbreak_individual
        add_outbreak_individual(camp, start_day=1, target_num_individuals=IMPORT_PER_MONTH, node_ids=[NODE],
                                repetitions=SIM_DAYS // 30, timesteps_between_repetitions=30, ignore_immunity=False)
        print(f"IMPORTATION: {IMPORT_PER_MONTH} infections every 30 days into node {NODE}")
    fix_property_restrictions(camp.campaign_dict)
    return camp


def update_serialize_parameters(simulation, df, x: int):
    path = df["serialized_file_path"][x]
    xTLH = float(df["x_Temporary_Larval_Habitat"][x])
    seed = int(df["Run_Number"][x])
    simulation.task.set_parameter("Serialized_Population_Filenames", df["Serialized_Population_Filenames"][x])
    simulation.task.set_parameter("Serialized_Population_Path", os.path.join(path, "output"))
    simulation.task.set_parameter("x_Temporary_Larval_Habitat", xTLH)
    simulation.task.set_parameter("Run_Number", seed)
    return {"xTLH": xTLH, "Run_Number": seed}


def build_demog():
    nodes = pd.read_csv(NODES_FILE)
    nodes["pop"] = TEST_POP
    nodes.to_csv("lake_nodes_bungoma_test.csv", index=False)
    demog = Demographics.from_csv(input_file="lake_nodes_bungoma_test.csv", id_ref=ID_REF,
                                  init_prev=0.01, include_biting_heterogeneity=True)
    demog.SetEquilibriumVitalDynamics()
    demog.SetAgeDistribution(Distributions.AgeDistribution_SSAfrica)
    return demog


def general_sim(selected_platform):
    platform = Platform(selected_platform, job_directory=manifest.job_directory, partition=manifest.partition,
                        time="12:00:00", modules=[manifest.singularity_module],
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
        raise SystemExit("Missing intervention files:\n  " + "\n  ".join(miss))
    parts = [build_burnin_df(i, platform, SERIALIZE_YEARS * 365) for i in BURNIN_EXP_IDS]
    burnin_df = pd.concat(parts, ignore_index=True)
    has = [all(os.path.exists(os.path.join(burnin_df["serialized_file_path"][i], "output", fn))
               for fn in burnin_df["Serialized_Population_Filenames"][i]) for i in burnin_df.index]
    burnin_df = burnin_df[has].copy()
    burnin_df["x_Temporary_Larval_Habitat"] = burnin_df["x_Temporary_Larval_Habitat"].astype(float).round(4)
    burnin_df = burnin_df.drop_duplicates(subset=["x_Temporary_Larval_Habitat", "Run_Number"]).reset_index(drop=True)
    if len(burnin_df) != EXPECTED_STATES:
        raise SystemExit(f"Expected {EXPECTED_STATES} burn-in states, found {len(burnin_df)}. Stopping.")
    if SUBSET:
        burnin_df = burnin_df[burnin_df["x_Temporary_Larval_Habitat"].isin([8.9797, 16.0]) & burnin_df["Run_Number"].isin([0, 1])].reset_index(drop=True)
    print("SCENARIO:", SCENARIO, "| IRS_REAL:", IRS_REAL, "| IMPORT per 30 days:", IMPORT_PER_MONTH, "| subset:", SUBSET)
    print("Burn-in simulations found:", len(burnin_df))
    print(burnin_df[["simid", "Run_Number", "x_Temporary_Larval_Habitat", "Serialized_Population_Filenames"]].head())
    builder = SimulationBuilder()
    builder.add_sweep_definition(partial(update_serialize_parameters, df=burnin_df), range(len(burnin_df.index)))
    end_day = SIM_DAYS
    add_report_event_counter(task, manifest, start_day=1, end_day=end_day, node_ids=[NODE],
                             min_age_years=0, max_age_years=100, event_trigger_list=EVENTS_TO_RECORD,
                             filename_suffix="node_" + str(NODE))
    add_malaria_summary_report(task, manifest, start_day=1, end_day=end_day, reporting_interval=30,
                               age_bins=[0.25, 5, 115], max_number_reports=SIM_DAYS // 30 + 2,
                               filename_suffix="monthly", pretty_format=True)
    add_spatial_report_malaria_filtered(task, manifest, start_day=1, end_day=end_day, reporting_interval=1,
                                        node_ids=None, min_age_years=0.25, max_age_years=100,
                                        spatial_output_channels=CHANNELS, filename_suffix="all_ages")
    add_spatial_report_malaria_filtered(task, manifest, start_day=1, end_day=end_day, reporting_interval=1,
                                        node_ids=None, min_age_years=0.5, max_age_years=5,
                                        spatial_output_channels=CHANNELS, filename_suffix="children")
    add_report_vector_stats(task, manifest, species_list=SPECIES, stratify_by_species=True)
    experiment = Experiment.from_builder(builder, task, name=f"{os.getlogin()}_LAKE_BUNGOMA_pickup_v3_nasa_{SCENARIO}_imp{IMPORT_PER_MONTH}_net{NET_TAG}")
    experiment.run(wait_until_done=False, platform=platform)
    print("\n=== BUNGOMA PICKUP SUBMITTED ===")
    print("Pickup ID:", experiment.uid)


if __name__ == "__main__":
    import emod_malaria.bootstrap as dtk
    dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    os.chmod(manifest.eradication_path, 0o755)
    general_sim("SLURM_LOCAL")
