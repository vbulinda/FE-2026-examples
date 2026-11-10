"""
Stage 3: run_inputs.py  (Lake Victoria, 8 counties)
Adapted from your working run_example.py plus the tutorial's Week 2 inputs script
and the Week 4 spatial build_demog(). One year, one seed, small population.

Expects (under manifest.input_dir):
  demographics/lake_nodes.csv                              node_id, lat, lon, pop  (ids 1-8)
  climate/LAKE/2010001-2019365/dtk_15arcmin_*_daily.bin    + .json (IdReference = LAKE_VICTORIA)
  vectors.csv   (optional; from fetch_vectors_lake.py)     used only for the species list
"""
import os
import pathlib
import pandas as pd

from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment

from emodpy.emod_task import EMODTask
import emod_api.campaign as camp

import emodpy_malaria.demographics.MalariaDemographics as Demographics
import emod_api.demographics.PreDefinedDistributions as Distributions

import manifest

# ---- settings for this stage -------------------------------------------------
ID_REF       = "LAKE_VICTORIA"            # must equal IdReference inside the climate .json files
SIM_YEARS    = 1
RUN_NUMBER   = 1
TEST_POP     = 1000                        # people per node for this test (real runs use lake_nodes.csv pop)
CLIMATE_DIR  = os.path.join("climate", "LAKE", "2010001-2019365")   # path inside the simulation
DEFAULT_SPECIES = ["gambiae", "arabiensis", "funestus"]
NODES_FILE   = os.path.join(manifest.input_dir, "demographics", "lake_nodes.csv")

# Centroids of the 8 Lake-endemic counties (same ids as fetch_climate_lake.py)
COUNTY_NODES = {
    "Kisumu":   {"id": 1, "lat": -0.1022, "lon": 34.7617},
    "Siaya":    {"id": 2, "lat":  0.0611, "lon": 34.2882},
    "Migori":   {"id": 3, "lat": -0.9995, "lon": 34.4744},
    "Homa Bay": {"id": 4, "lat": -0.5284, "lon": 34.4571},
    "Kakamega": {"id": 5, "lat":  0.2842, "lon": 34.7533},
    "Busia":    {"id": 6, "lat":  0.4601, "lon": 34.1115},
    "Bungoma":  {"id": 7, "lat":  0.5695, "lon": 34.5606},
    "Vihiga":   {"id": 8, "lat":  0.0815, "lon": 34.7224},
}
CLIMATE_FILES = ["dtk_15arcmin_air_temperature_daily", "dtk_15arcmin_land_temperature_daily",
                 "dtk_15arcmin_rainfall_daily", "dtk_15arcmin_relative_humidity_daily"]


def preflight():
    """Create the node table if missing; stop early with a clear message if climate files are missing."""
    if not os.path.exists(NODES_FILE):
        os.makedirs(os.path.dirname(NODES_FILE), exist_ok=True)
        rows = [{"node_id": v["id"], "lat": v["lat"], "lon": v["lon"], "pop": 10000} for v in COUNTY_NODES.values()]
        pd.DataFrame(rows).sort_values("node_id").to_csv(NODES_FILE, index=False)
        print("Wrote", NODES_FILE)
    cdir = os.path.join(manifest.input_dir, CLIMATE_DIR)
    missing = [f for f in CLIMATE_FILES for ext in (".bin", ".bin.json")
               if not os.path.exists(os.path.join(cdir, f + ext))]
    if missing:
        raise SystemExit("Climate files missing in %s:\n  %s\nRun: python fetch_climate_lake.py" % (cdir, "\n  ".join(missing)))
    print("Preflight OK: 8 nodes, climate files found in", cdir)
VECTOR_FILE  = os.path.join(manifest.input_dir, "vectors.csv")


def species_list():
    """Species from vectors.csv if present, otherwise the three default Anopheles."""
    if os.path.exists(VECTOR_FILE):
        v = pd.read_csv(VECTOR_FILE)
        s = [x for x in v["species"].unique()]
        if s:
            return s
    return DEFAULT_SPECIES


def set_param_fn(config):
    import emodpy_malaria.malaria_config as conf
    config = conf.set_team_defaults(config, manifest)

    # Needed so InsetChart.json is written (off by default in the installed build)
    config.parameters.Enable_Default_Reporting = 1

    config.parameters.Simulation_Duration = SIM_YEARS * 365
    config.parameters.Run_Number = RUN_NUMBER

    # Vectors (species are global in EMOD; node-level differences come later via habitat scaling)
    conf.add_species(config, manifest, species_list())

    # Climate from your fetched files
    config.parameters.Climate_Model = "CLIMATE_BY_DATA"
    config.parameters.Air_Temperature_Filename      = os.path.join(CLIMATE_DIR, "dtk_15arcmin_air_temperature_daily.bin")
    config.parameters.Land_Temperature_Filename     = os.path.join(CLIMATE_DIR, "dtk_15arcmin_land_temperature_daily.bin")
    config.parameters.Rainfall_Filename             = os.path.join(CLIMATE_DIR, "dtk_15arcmin_rainfall_daily.bin")
    config.parameters.Relative_Humidity_Filename    = os.path.join(CLIMATE_DIR, "dtk_15arcmin_relative_humidity_daily.bin")
    return config


def build_camp():
    camp.set_schema(manifest.schema_file)
    return camp


def build_demog():
    """8 nodes from lake_nodes.csv, population scaled down to TEST_POP for the test."""
    nodes = pd.read_csv(NODES_FILE)
    assert list(nodes.columns[:4]) == ["node_id", "lat", "lon", "pop"], "lake_nodes.csv needs node_id,lat,lon,pop"
    nodes["pop"] = TEST_POP
    test_csv = "lake_nodes_test.csv"
    nodes.to_csv(test_csv, index=False)

    demog = Demographics.from_csv(input_file=test_csv, id_ref=ID_REF,
                                  init_prev=0.01, include_biting_heterogeneity=True)
    demog.SetEquilibriumVitalDynamics()
    demog.SetAgeDistribution(Distributions.AgeDistribution_SSAfrica)
    return demog


def general_sim(selected_platform):
    platform = Platform(selected_platform, job_directory=manifest.job_directory,
                        partition=manifest.partition, time=manifest.sim_time,
                        modules=[manifest.singularity_module],
                        max_running_jobs=manifest.max_running_jobs)

    print("Creating EMODTask (from files)...")
    task = EMODTask.from_default2(
        config_path="config.json",
        eradication_path=manifest.eradication_path,
        campaign_builder=build_camp,
        schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn,
        ep4_custom_cb=None,
        demog_builder=build_demog,
        plugin_report=None,
    )
    task.set_sif(manifest.SIF_PATH, platform)

    # climate folder goes in as an asset, same relative path the config uses
    task.common_assets.add_directory(os.path.join(manifest.input_dir, "climate"), relative_path="climate")

    experiment = Experiment.from_task(task, name="lake_inputs_test")
    experiment.run(wait_until_done=True, platform=platform)

    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.\n")
        exit()
    print(f"Experiment {experiment.uid} succeeded.")


if __name__ == "__main__":
    import emod_malaria.bootstrap as dtk
    preflight()
    dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    os.chmod(manifest.eradication_path, 0o755)

    general_sim("SLURM_LOCAL")