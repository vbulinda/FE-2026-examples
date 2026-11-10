"""
run_bungoma_burnin_v2.py - Bungoma (node 7) burn-in v2, 50 years, WITH baseline case
management (the 2010 rows of case_management_v2.csv) and importation from day 1,
as in the Migori and Homa Bay burn-ins. Built from run_homabay_burnin_v2.py.
Species: Bungoma values of burn-in v1 (gambiae x20, arabiensis x0.0225, funestus x0.006 -> 97.2/2.2/0.6 %).
xTLH grid 0.5 to 8 (13 values, geomspace) x 4 seeds = 52 runs (v1 grid 1 to 32 had its best value near the low edge).
Climate: the 10-year 2010-2019 file LAKE_KISUMU (3,650 days; 18,250 = 5 x 3,650 keeps the season aligned).
Run from FE-2026-examples:   python run_bungoma_burnin_v2.py     It prints the experiment ID and stops.
"""
import os, pathlib
import numpy as np
import pandas as pd
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from emodpy.emod_task import EMODTask
import emod_api.campaign as camp
import emodpy_malaria.demographics.MalariaDemographics as Demographics
import emod_api.demographics.PreDefinedDistributions as Distributions
from emodpy_malaria.reporters.builtin import *
import emodpy_malaria.interventions.treatment_seeking as cm
import manifest

ID_REF = "LAKE_VICTORIA"
SERIALIZE_YEARS = 50
NODE = 1
IMPORT_PER_MONTH = int(os.environ.get("IMPORT", "1"))
INTERV_DIR = os.path.join(manifest.input_dir, "interventions")
sim_years = SERIALIZE_YEARS
SEEDS = [int(v) for v in os.environ.get("SEEDS", "0,1,2,3").split(",")]
SUF = os.environ.get("SUF", "")
XTLH_VALUES = [float(v) for v in os.environ["XTLH_LIST"].split(",")] if os.environ.get("XTLH_LIST") else [round(float(x), 4) for x in np.geomspace(1.0, 16.0, 13)]
TEST_POP = 5000                       # same as the earlier calibration burn-ins; the pickup must use the same value
RUN_TIME = "12:00:00"
CLIMATE_DIR = os.path.join("climate", "LAKE_KISUMU", "2010001-2019365")
NODES_FILE = os.path.join(manifest.input_dir, "demographics", "lake_nodes_kisumu.csv")
SPECIES = ["gambiae", "arabiensis", "funestus"]
MIX = os.environ.get("MIX", "K")
G_MULT = float(os.environ.get("G_MULT", "7.527"))
ARA_MULT = float(os.environ.get("ARA_MULT", "0.5751"))
FUN_MULT = float(os.environ.get("FUN_MULT", "0.2235"))   # Bungoma species test, log step 43e
FILES = {"air": "dtk_15arcmin_air_temperature_daily", "land": "dtk_15arcmin_land_temperature_daily",
         "rh": "dtk_15arcmin_relative_humidity_daily", "rain": "dtk_15arcmin_rainfall_daily"}


def check_inputs():
    cdir = os.path.join(manifest.input_dir, CLIMATE_DIR)
    need = [NODES_FILE] + [os.path.join(cdir, n + e) for n in FILES.values() for e in (".bin", ".bin.json")]
    miss = [p for p in need if not os.path.exists(p)]
    if miss:
        raise SystemExit("Missing input files (run python make_kisumu_inputs.py first):\n  " + "\n  ".join(miss))


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
    conf.set_max_larval_capacity(config, "gambiae", "WATER_VEGETATION", int(20000000 * G_MULT))
    conf.set_max_larval_capacity(config, "arabiensis", "TEMPORARY_RAINFALL", int(800000000 * ARA_MULT))
    conf.set_max_larval_capacity(config, "arabiensis", "CONSTANT", int(80000000 * ARA_MULT))
    conf.set_max_larval_capacity(config, "funestus", "TEMPORARY_RAINFALL", int(800000000 * FUN_MULT))
    conf.set_max_larval_capacity(config, "funestus", "CONSTANT", int(80000000 * FUN_MULT))
    config.parameters.Simulation_Duration = sim_years * 365
    config.parameters.Run_Number = 0
    config.parameters.x_Temporary_Larval_Habitat = 1.0
    config.parameters.Serialized_Population_Writing_Type = "TIMESTEP"
    config.parameters.Serialization_Time_Steps = [365 * SERIALIZE_YEARS]
    config.parameters.Serialization_Mask_Node_Write = 0
    config.parameters.Serialization_Precision = "REDUCED"
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


def build_camp():
    camp.set_schema(manifest.schema_file)
    cm_all = pd.read_csv(os.path.join(INTERV_DIR, "case_management_v2.csv"))
    sel = (cm_all["node_id"].astype(int) == NODE) & (cm_all["phase"] == "pickup")
    sub = cm_all[sel & (cm_all["year"] == 2010)].reset_index(drop=True)
    assert len(sub) == 4, "expected 4 case-management rows for 2010"
    targets = [{"trigger": str(sub["trigger"][r]), "coverage": float(sub["coverage"][r]),
                "agemin": float(sub["age_min"][r]), "agemax": float(sub["age_max"][r]),
                "rate": float(sub["rate"][r])} for r in range(len(sub))]
    cm.add_treatment_seeking(camp, node_ids=[NODE], start_day=1,
                             duration=SERIALIZE_YEARS * 365, drug=["Artemether", "Lumefantrine"],
                             targets=targets, broadcast_event_name="Received_Treatment")
    print("BASELINE case management (2010 rows) applied for", SERIALIZE_YEARS, "years")
    if IMPORT_PER_MONTH > 0:
        from emodpy_malaria.interventions.outbreak import add_outbreak_individual
        add_outbreak_individual(camp, start_day=1, target_num_individuals=IMPORT_PER_MONTH,
                                node_ids=[NODE], repetitions=SERIALIZE_YEARS * 365 // 30,
                                timesteps_between_repetitions=30, ignore_immunity=False)
        print("IMPORTATION:", IMPORT_PER_MONTH, "infections every 30 days")
    fix_property_restrictions(camp.campaign_dict)
    return camp


def build_demog():
    nodes = pd.read_csv(NODES_FILE)
    assert list(nodes.columns[:4]) == ["node_id", "lat", "lon", "pop"]
    nodes["pop"] = TEST_POP
    nodes.to_csv("lake_nodes_kisumu_test.csv", index=False)      # own file name: does not overwrite the 8-node test file
    demog = Demographics.from_csv(input_file="lake_nodes_kisumu_test.csv", id_ref=ID_REF,
                                  init_prev=0.01, include_biting_heterogeneity=True)
    demog.SetEquilibriumVitalDynamics()
    demog.SetAgeDistribution(Distributions.AgeDistribution_SSAfrica)
    return demog


def set_pair(simulation, pair):
    x, s = pair
    simulation.task.set_parameter("x_Temporary_Larval_Habitat", x)
    simulation.task.set_parameter("Run_Number", s)
    return {"x_Temporary_Larval_Habitat": x, "Run_Number": s}


def general_sim(selected_platform):
    platform = Platform(selected_platform, job_directory=manifest.job_directory, partition=manifest.partition,
                        time=RUN_TIME, modules=[manifest.singularity_module], max_running_jobs=manifest.max_running_jobs)
    task = EMODTask.from_default2(config_path="config.json", eradication_path=manifest.eradication_path,
                                  campaign_builder=build_camp, schema_path=manifest.schema_file,
                                  param_custom_cb=set_param_fn, ep4_custom_cb=None,
                                  demog_builder=build_demog, plugin_report=None)
    task.set_sif(manifest.SIF_PATH, platform)
    task.common_assets.add_directory(os.path.join(manifest.input_dir, "climate"), relative_path="climate")
    end_day = sim_years * 365
    add_malaria_summary_report(task, manifest, start_day=1, end_day=end_day, reporting_interval=30,
                               age_bins=[0.25, 5, 115], max_number_reports=SERIALIZE_YEARS * 13,
                               filename_suffix="monthly", pretty_format=True)
    chans = ["Population", "PCR_Parasite_Prevalence", "Blood_Smear_Parasite_Prevalence", "New_Clinical_Cases"]
    add_spatial_report_malaria_filtered(task, manifest, start_day=1, end_day=end_day, reporting_interval=1,
                                        node_ids=None, min_age_years=0.25, max_age_years=100,
                                        spatial_output_channels=chans, filename_suffix="all_ages")
    # per-species vector counts every day: measures the species mix the model produces
    add_report_vector_stats(task, manifest, species_list=SPECIES, stratify_by_species=True)
    user = os.getlogin()
    builder = SimulationBuilder()
    pairs = [(x, s) for x in XTLH_VALUES for s in SEEDS]
    builder.add_sweep_definition(set_pair, pairs)
    experiment = Experiment.from_builder(builder, task, name=f"{user}_LAKE_KISUMU_burnin_{MIX}{SUF}")
    experiment.run(wait_until_done=False, platform=platform)
    print("\n=== KISUMU BURN-IN v2 SUBMITTED ===")
    print("Mix", MIX, "G", G_MULT, "ARA", ARA_MULT, "FUN", FUN_MULT)
    print("Burn-in ID:", experiment.uid, " simulations:", len(pairs))
    print("xTLH values:", XTLH_VALUES)
    print("Write this ID in the log. Check with: idmtools slurm", manifest.job_directory, "status --exp-id", experiment.uid)


if __name__ == "__main__":
    import emod_malaria.bootstrap as dtk
    check_inputs()
    dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    os.chmod(manifest.eradication_path, 0o755)
    general_sim("SLURM_LOCAL")