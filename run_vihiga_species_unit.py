"""
run_vihiga_species_unit.py - vector-only test, Vihiga node, 3 years, one experiment per multiplier set (2 experiments of one simulation).
Question: which habitat multipliers (gambiae WATER_VEGETATION; arabiensis and funestus TEMPORARY_RAINFALL and CONSTANT)
reproduce (A) a funestus-dominant mix (about 85%, as in Migori before spraying, Abong'o 2020) and (B) an arabiensis/gambiae-dominant
mix (about 49% arabiensis, 32% gambiae s.s., 19% funestus, as at Kisian, Kisumu County, 2018) at Homa Bay climate?
Grid: G_MULT in (1, 4, 16) x ARA_MULT in (1, 0.25, 0.0625) x FUN_MULT in (1, 0.25, 0.0625), xTLH 4.0, no interventions.
No saved state is written (this is not a burn-in). Writes vihiga_species_ids.txt (g, ara, fun, experiment id).
Run from FE-2026-examples:   python run_vihiga_species_unit.py
"""
import os, pathlib
import pandas as pd
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from emodpy.emod_task import EMODTask
import emod_api.campaign as camp
import emodpy_malaria.demographics.MalariaDemographics as Demographics
import emod_api.demographics.PreDefinedDistributions as Distributions
from emodpy_malaria.reporters.builtin import *
import manifest

ID_REF = "LAKE_VICTORIA"
SIM_YEARS = 3
XTLH = 4.0
GRID = [(1, 1.0, 1.0), (1, 0.125, 1.0)]
TEST_POP = 5000
RUN_TIME = "04:00:00"
CLIMATE_DIR = os.path.join("climate", "LAKE_VIHIGA", "2010001-2019365")
NODES_FILE = os.path.join(manifest.input_dir, "demographics", "lake_nodes_vihiga.csv")
SPECIES = ["gambiae", "arabiensis", "funestus"]
FILES = {"air": "dtk_15arcmin_air_temperature_daily", "land": "dtk_15arcmin_land_temperature_daily",
         "rh": "dtk_15arcmin_relative_humidity_daily", "rain": "dtk_15arcmin_rainfall_daily"}


def make_set_param_fn(trio):
    g_mult, ara_mult, fun_mult = trio
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
        conf.set_max_larval_capacity(config, "gambiae", "WATER_VEGETATION", int(20000000 * g_mult))
        conf.set_max_larval_capacity(config, "arabiensis", "TEMPORARY_RAINFALL", int(800000000 * ara_mult))
        conf.set_max_larval_capacity(config, "arabiensis", "CONSTANT", int(80000000 * ara_mult))
        conf.set_max_larval_capacity(config, "funestus", "TEMPORARY_RAINFALL", int(800000000 * fun_mult))
        conf.set_max_larval_capacity(config, "funestus", "CONSTANT", int(80000000 * fun_mult))
        config.parameters.Simulation_Duration = SIM_YEARS * 365
        config.parameters.Run_Number = 0
        config.parameters.x_Temporary_Larval_Habitat = XTLH
        return config
    return set_param_fn


def build_camp():
    camp.set_schema(manifest.schema_file)
    return camp


def build_demog():
    nodes = pd.read_csv(NODES_FILE)
    nodes["pop"] = TEST_POP
    nodes.to_csv("lake_nodes_vihiga_test.csv", index=False)
    demog = Demographics.from_csv(input_file="lake_nodes_vihiga_test.csv", id_ref=ID_REF,
                                  init_prev=0.01, include_biting_heterogeneity=True)
    demog.SetEquilibriumVitalDynamics()
    demog.SetAgeDistribution(Distributions.AgeDistribution_SSAfrica)
    return demog


def tag_fn(simulation, v):
    return {"g_mult": v[0], "ara_mult": v[1], "fun_mult": v[2]}


def run_one(trio, platform):
    task = EMODTask.from_default2(config_path="config.json", eradication_path=manifest.eradication_path,
                                  campaign_builder=build_camp, schema_path=manifest.schema_file,
                                  param_custom_cb=make_set_param_fn(trio), ep4_custom_cb=None,
                                  demog_builder=build_demog, plugin_report=None)
    task.set_sif(manifest.SIF_PATH, platform)
    task.common_assets.add_directory(os.path.join(manifest.input_dir, "climate"), relative_path="climate")
    add_report_vector_stats(task, manifest, species_list=SPECIES, stratify_by_species=True)
    builder = SimulationBuilder()
    builder.add_sweep_definition(tag_fn, [trio])
    experiment = Experiment.from_builder(builder, task, name=f"{os.getlogin()}_LAKE_VIHIGA_species_g{trio[0]}_a{trio[1]}_f{trio[2]}")
    experiment.run(wait_until_done=False, platform=platform)
    print(f"G {trio[0]} ARA {trio[1]} FUN {trio[2]}: experiment {experiment.uid}")
    return str(experiment.uid)


if __name__ == "__main__":
    import emod_malaria.bootstrap as dtk
    dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    os.chmod(manifest.eradication_path, 0o755)
    platform = Platform("SLURM_LOCAL", job_directory=manifest.job_directory, partition=manifest.partition,
                        time=RUN_TIME, modules=[manifest.singularity_module], max_running_jobs=manifest.max_running_jobs)
    ids = {m: run_one(m, platform) for m in GRID}
    with open("vihiga_species_ids.txt", "w") as fh:
        for m, i in ids.items():
            fh.write(f"{m[0]},{m[1]},{m[2]},{i}\n")
    print("\n=== VIHIGA SPECIES GRID SUBMITTED: 2 experiments, ids in vihiga_species_ids.txt ===")
