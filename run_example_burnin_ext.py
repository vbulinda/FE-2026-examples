"""
run_lake_burnin_cal_rerun.py  -  reruns ONLY the calibration burn-in runs that have no saved state.

It looks at experiment 10b92b9f... , lists every (habitat value, seed) that has state-10950.dtk, compares with the
20 planned runs, and submits only the missing ones (a new experiment). Print the new ID and stop.
Run from FE-2026-examples:   python run_lake_burnin_cal_rerun.py

(original docstring follows)
run_lake_burnin30.py  -  ONE self-contained script (Stage 7: REAL burn-in, Lake Victoria, 8 counties)

What it does:
  Runs a 30-year burn-in (the population builds up immunity and reaches a stable starting state)
  for 12 values of x_Temporary_Larval_Habitat (0.001 to about 3) and 2 seeds = 24 simulations,
  then saves the population state at the end of year 30 for the pickup stage.

Why: the 2-year test burn-in gave 59-88% prevalence and 2,300-3,500 cases per 1,000 per year.
Two years is too short for immunity to develop, and 0.32 was already too high, so the range
now goes much lower. Longer burn-in plus lower habitat is expected to bring these down.

The script does NOT wait for the runs: it prints the burn-in ID and stops. Check with
  idmtools slurm <job_directory> status --exp-id <ID>
If runs hit the time limit, raise manifest.sim_time (for example '12:00:00').
Climate: the 2010-2019 data are expected to repeat for the 30 years (check Assets/ output if worried).

Run from the FE-2026-examples folder:   python run_lake_burnin30.py
"""
import pathlib
import os
from functools import partial
import numpy as np
import pandas as pd
from idmtools.builders import SimulationBuilder

#idmtools
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment

#emodpy
from emodpy.emod_task import EMODTask
import emod_api.campaign as camp

#emodpy-malaria
import emodpy_malaria.demographics.MalariaDemographics as Demographics
import emod_api.demographics.PreDefinedDistributions as Distributions
from emodpy_malaria.reporters.builtin import *

import manifest

# ======================= PART 1: climate builder (from fetch_climate_lake.py) =======================
import argparse, json, os, sys, time, datetime as dt
import numpy as np
import pandas as pd

ID_REF = "LAKE_VICTORIA"
POP = 10000                      # simulated people per node (use 1000 for quick tests; change in lake_nodes.csv)
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
FILES = {"air": "dtk_15arcmin_air_temperature_daily", "land": "dtk_15arcmin_land_temperature_daily",
         "rh": "dtk_15arcmin_relative_humidity_daily", "rain": "dtk_15arcmin_rainfall_daily"}

# ------------------------------------------------------------------ fetching
def _get(url, params, tries=4):
    import requests
    for k in range(tries):
        try:
            r = requests.get(url, params=params, timeout=90)
            if r.status_code == 200:
                return r.json()
            print(f"  HTTP {r.status_code}: {r.text[:150]}")
        except Exception as e:
            print(f"  request error: {e}")
        time.sleep(5 * (k + 1))
    sys.exit("Could not reach the data service. Check internet access from this machine, or try --source openmeteo.")

def fetch_power(lat, lon, y0, y1):
    frames = []
    for y in range(y0, y1 + 1):
        j = _get("https://power.larc.nasa.gov/api/temporal/daily/point",
                 {"parameters": "T2M,RH2M,PRECTOTCORR", "community": "AG", "longitude": lon, "latitude": lat,
                  "start": f"{y}0101", "end": f"{y}1231", "format": "JSON"})
        p = j["properties"]["parameter"]
        df = pd.DataFrame({"date": pd.to_datetime(list(p["T2M"].keys()), format="%Y%m%d"),
                           "temp": list(p["T2M"].values()), "rh": list(p["RH2M"].values()),
                           "rain": list(p["PRECTOTCORR"].values())})
        frames.append(df)
    out = pd.concat(frames).replace(-999, np.nan)
    out["rh"] = out["rh"] / 100.0
    return out

def fetch_openmeteo(lat, lon, y0, y1):
    j = _get("https://archive-api.open-meteo.com/v1/archive",
             {"latitude": lat, "longitude": lon, "start_date": f"{y0}-01-01", "end_date": f"{y1}-12-31",
              "daily": "temperature_2m_mean,relative_humidity_2m_mean,precipitation_sum", "timezone": "UTC"})
    d = j["daily"]
    return pd.DataFrame({"date": pd.to_datetime(d["time"]), "temp": d["temperature_2m_mean"],
                         "rh": np.array(d["relative_humidity_2m_mean"], dtype=float) / 100.0,
                         "rain": d["precipitation_sum"]})

def synthetic(lat, lon, y0, y1, seed):
    rng = np.random.default_rng(seed)
    dates = pd.date_range(f"{y0}-01-01", f"{y1}-12-31")
    doy = dates.dayofyear.values
    temp = 22 + 2.5 * np.sin(2 * np.pi * (doy - 60) / 365) + rng.normal(0, .6, len(doy))
    rain = np.clip(rng.gamma(0.5, 8, len(doy)) * (1 + np.sin(2 * np.pi * (doy - 90) / 182)), 0, None)
    rh = np.clip(0.75 + 0.1 * np.sin(2 * np.pi * (doy - 100) / 182) + rng.normal(0, .03, len(doy)), .3, .99)
    return pd.DataFrame({"date": dates, "temp": temp, "rh": rh, "rain": rain})

# ------------------------------------------------------------------ cleaning and writing
def to_365(df, y0, y1):
    df = df.drop_duplicates("date").set_index("date").sort_index()
    full = pd.date_range(f"{y0}-01-01", f"{y1}-12-31")
    df = df.reindex(full)
    miss = df.isna().mean().max()
    if miss > 0.02:
        sys.exit(f"{miss:.1%} of days are missing: too many to interpolate. Try the other --source.")
    df = df.interpolate(limit=10, limit_direction="both")
    df = df[~((df.index.month == 2) & (df.index.day == 29))]
    df["rain"] = df["rain"].clip(lower=0)
    df["rh"] = df["rh"].clip(0.0, 1.0)
    assert len(df) == 365 * (y1 - y0 + 1), len(df)
    return df

def write_bin(outdir, stem, arr, node_ids, y0, y1, source, lat, lon):
    """arr shape (n_nodes, n_days), written node-major as little-endian float32 (EMOD schema 2.0)."""
    n, days = arr.shape
    arr.astype("<f4").tofile(os.path.join(outdir, stem + ".bin"))
    meta = {"Metadata": {
        "Tool": "fetch_climate_lake.py", "Author": "Lake Victoria EMOD project", "WeatherSchemaVersion": "2.0",
        "IdReference": ID_REF, "WeatherCellCount": n, "NodeCount": n, "OffsetEntryCount": n,
        "NumberDTKNodes": n, "DatavaluePerCell": days, "DatavalueCount": days,
        "UpdateResolution": "CLIMATE_UPDATE_DAY", "OriginalDataYears": f"{y0}-{y1}",
        "StartDayOfYear": "January 1", "DataProvenance": source, "Resolution": "point (county centroid)",
        "UpperLatitude": max(lat), "LeftLongitude": min(lon), "BottomLatitude": min(lat), "RightLongitude": max(lon)},
        "NodeOffsets": "".join(f"{nid:08x}{i * days * 4:08x}" for i, nid in enumerate(node_ids))}
    with open(os.path.join(outdir, stem + ".bin.json"), "w") as f:
        json.dump(meta, f, indent=4)

def read_bin(outdir, stem):
    with open(os.path.join(outdir, stem + ".bin.json")) as f:
        m = json.load(f)
    days, n = m["Metadata"]["DatavalueCount"], m["Metadata"]["NodeCount"]
    ids = [int(m["NodeOffsets"][16 * i:16 * i + 8], 16) for i in range(n)]
    offs = [int(m["NodeOffsets"][16 * i + 8:16 * i + 16], 16) for i in range(n)]
    raw = np.fromfile(os.path.join(outdir, stem + ".bin"), dtype="<f4")
    assert raw.size == n * days, f"{stem}: expected {n*days} values, found {raw.size}"
    return ids, np.array([raw[o // 4:o // 4 + days] for o in offs]), m

def verify(outdir, plot=True):
    id2name = {c["id"]: n for n, c in COUNTY_NODES.items()}
    print(f"\nVerification of {outdir}")
    series = {}
    for key, stem in FILES.items():
        ids, a, m = read_bin(outdir, stem)
        series[key] = a
        print(f"  {stem}: {len(ids)} nodes x {a.shape[1]} days, IdReference={m['Metadata']['IdReference']}, "
              f"range {a.min():.2f} to {a.max():.2f}")
    print("\n  node  county     mean temp C   mean RH   rain mm/year")
    yrs = series["air"].shape[1] / 365
    for i, nid in enumerate(ids):
        print(f"  {nid:<5} {id2name.get(nid,'?'):<10} {series['air'][i].mean():>10.1f} {series['rh'][i].mean():>10.2f} "
              f"{series['rain'][i].sum() / yrs:>12.0f}")
    print("\n  Sanity: Lake Victoria basin should show about 20-25 C, RH 0.6-0.85, rain roughly 1000-2000 mm/year, "
          "with two rainy seasons (March-May and Oct-Dec).")
    if plot:
        try:
            import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
            fig, ax = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
            x = np.arange(365)
            for lab, k, a in (("Air temperature (C)", "air", ax[0]), ("Relative humidity", "rh", ax[1]),
                              ("Rainfall (mm/day)", "rain", ax[2])):
                mean_year = series[k].reshape(len(ids), -1, 365).mean(axis=1)
                for i, nid in enumerate(ids):
                    a.plot(x, pd.Series(mean_year[i]).rolling(15, center=True, min_periods=1).mean(), label=id2name.get(nid))
                a.set_ylabel(lab)
            ax[2].set_xlabel("Day of year (mean over all years, 15-day smooth)"); ax[0].legend(ncol=4, fontsize=7)
            fig.tight_layout(); p = os.path.join(outdir, "climate_check.png"); fig.savefig(p, dpi=120); print(f"  Plot: {p}")
        except Exception as e:
            print("  (plot skipped:", e, ")")

# ------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=2010); ap.add_argument("--end", type=int, default=2019)
    ap.add_argument("--source", choices=["power", "openmeteo"], default="power")
    ap.add_argument("--inputs", default="inputs", help="EMOD inputs folder (default: ./inputs)")
    ap.add_argument("--selftest", action="store_true"); ap.add_argument("--verify", action="store_true")
    a = ap.parse_args(argv)
    tag = f"{a.start}001-{a.end}365"
    outdir = os.path.join(a.inputs, "climate", "LAKE", tag)
    demo_dir = os.path.join(a.inputs, "demographics")
    if a.selftest:
        outdir = os.path.join("selftest_out", tag); demo_dir = os.path.join("selftest_out", "demographics")
    os.makedirs(outdir, exist_ok=True); os.makedirs(demo_dir, exist_ok=True)
    if a.verify:
        verify(outdir); return
    nodes = sorted(COUNTY_NODES.items(), key=lambda kv: kv[1]["id"])
    ids = [v["id"] for _, v in nodes]; lats = [v["lat"] for _, v in nodes]; lons = [v["lon"] for _, v in nodes]
    cache = os.path.join(outdir, "raw_downloads"); os.makedirs(cache, exist_ok=True)
    data = {k: [] for k in ("temp", "rh", "rain")}
    for name, v in nodes:
        f = os.path.join(cache, f"{name.replace(' ', '_')}_{a.source if not a.selftest else 'synthetic'}.csv")
        if os.path.exists(f):
            df = pd.read_csv(f, parse_dates=["date"]); print(f"{name}: using cached {f}")
        else:
            print(f"{name}: fetching {a.start}-{a.end} ({'synthetic' if a.selftest else a.source}) ...")
            df = synthetic(v["lat"], v["lon"], a.start, a.end, v["id"]) if a.selftest else \
                 (fetch_power if a.source == "power" else fetch_openmeteo)(v["lat"], v["lon"], a.start, a.end)
            df.to_csv(f, index=False)
        df = to_365(df, a.start, a.end)
        for k in data: data[k].append(df[k].values)
    arrs = {k: np.array(v) for k, v in data.items()}
    src = "SYNTHETIC TEST DATA, NOT REAL CLIMATE" if a.selftest else \
          ("NASA POWER daily point data (T2M, RH2M, PRECTOTCORR)" if a.source == "power" else "Open-Meteo archive (ERA5)")
    write_bin(outdir, FILES["air"], arrs["temp"], ids, a.start, a.end, src, lats, lons)
    write_bin(outdir, FILES["land"], arrs["temp"], ids, a.start, a.end, src, lats, lons)
    write_bin(outdir, FILES["rh"], arrs["rh"], ids, a.start, a.end, src, lats, lons)
    write_bin(outdir, FILES["rain"], arrs["rain"], ids, a.start, a.end, src, lats, lons)
    pd.DataFrame({"node_id": ids, "lat": lats, "lon": lons, "pop": POP}).to_csv(os.path.join(demo_dir, "lake_nodes.csv"), index=False)
    pd.DataFrame({"node_id": ids, "county": [n for n, _ in nodes], "lat": lats, "lon": lons}).to_csv(
        os.path.join(demo_dir, "lake_nodes_county_names.csv"), index=False)
    print(f"\nWrote 4 climate file pairs to {outdir}\nWrote {os.path.join(demo_dir, 'lake_nodes.csv')}")
    verify(outdir)


# ======================= PART 2: the EMOD run =======================
# ---------------- Lake Victoria settings ----------------
SERIALIZE_YEARS = 30                # real burn-in
sim_years   = SERIALIZE_YEARS
NUM_SEEDS   = 2
XTLH_VALUES = [1.9373, 2.5022, 3.2317, 4.1738]
OLD_EXP_PREFIX = "NO_SUCH_FOLDER_"
RERUN_TIME  = "12:00:00"           # time limit per run; if the cluster refuses it, try 08:00:00
TEST_POP    = 5000                  # people per node (5x the earlier burn-in, fewer chance extinctions). The pickup must use the same value
CLIMATE_DIR = os.path.join('climate', 'LAKE', '2010001-2019365')      # path inside the simulation
CLIMATE_TAG = '2010001-2019365'
NODES_FILE  = os.path.join(manifest.input_dir, 'demographics', 'lake_nodes.csv')
SPECIES     = ["gambiae", "arabiensis", "funestus"]
# ID_REF ("LAKE_VICTORIA") is defined in part 1 and must equal IdReference in the climate .json files
# --------------------------------------------------------


def ensure_inputs():
    """Build lake_nodes.csv and the climate files if they are not there yet."""
    cdir = os.path.join(manifest.input_dir, 'climate', 'LAKE', CLIMATE_TAG)
    needed = [os.path.join(cdir, n + ext) for n in FILES.values() for ext in ('.bin', '.bin.json')]
    if os.path.exists(NODES_FILE) and all(os.path.exists(p) for p in needed):
        print("Inputs found: lake_nodes.csv + 8 climate files.")
        return
    print("Inputs missing, building them now (needs internet access from this server) ...")
    last = None
    for src in ("power", "openmeteo"):
        try:
            print("\n=== Downloading climate from", src, "===")
            main(["--source", src, "--inputs", manifest.input_dir])
            break
        except BaseException as e:          # fetch code exits with sys.exit on network failure
            last = e
            print("!!", src, "failed:", type(e).__name__, e)
    missing = [p for p in needed if not os.path.exists(p)]
    if missing or not os.path.exists(NODES_FILE):
        raise SystemExit("Could not build the inputs. Last error: %s\nPaste the messages above to Claude." % last)


def set_param_fn(config):
    import emodpy_malaria.malaria_config as conf
    config = conf.set_team_defaults(config, manifest)
    config.parameters.Climate_Model = "CLIMATE_BY_DATA"

    # Defaults to 0 (off) in the installed build - turn on so InsetChart.json is written
    config.parameters.Enable_Default_Reporting = 1

    config.parameters.Air_Temperature_Filename   = os.path.join(CLIMATE_DIR, FILES["air"] + '.bin')
    config.parameters.Land_Temperature_Filename  = os.path.join(CLIMATE_DIR, FILES["land"] + '.bin')
    config.parameters.Rainfall_Filename          = os.path.join(CLIMATE_DIR, FILES["rain"] + '.bin')
    config.parameters.Relative_Humidity_Filename = os.path.join(CLIMATE_DIR, FILES["rh"] + '.bin')

    conf.add_species(config, manifest, SPECIES)

    config.parameters.Simulation_Duration = sim_years*365
    config.parameters.Run_Number = 0

    # Larval habitat multiplier: the parameter that will be calibrated (overridden by the sweep below)
    config.parameters.x_Temporary_Larval_Habitat = 1.0

    # Serialization: write the population state at the end of the burn-in
    config.parameters.Serialized_Population_Writing_Type = "TIMESTEP"
    config.parameters.Serialization_Time_Steps = [365 * SERIALIZE_YEARS]
    config.parameters.Serialization_Mask_Node_Write = 0
    config.parameters.Serialization_Precision = "REDUCED"
    return config


def set_param(simulation, param, value):
    """Sets one config parameter in one simulation; used by the builder."""
    return simulation.task.set_parameter(param, value)


def build_camp():
    camp.set_schema(manifest.schema_file)
    return camp


def build_demog():
    """8 Lake Victoria counties from lake_nodes.csv (population reduced to TEST_POP for the test)."""
    nodes = pd.read_csv(NODES_FILE)
    assert list(nodes.columns[:4]) == ["node_id", "lat", "lon", "pop"], "lake_nodes.csv needs node_id,lat,lon,pop"
    nodes['pop'] = TEST_POP
    nodes.to_csv('lake_nodes_test.csv', index=False)
    demog = Demographics.from_csv(input_file='lake_nodes_test.csv', id_ref=ID_REF,
                                  init_prev=0.01, include_biting_heterogeneity=True)
    demog.SetEquilibriumVitalDynamics()
    demog.SetAgeDistribution(Distributions.AgeDistribution_SSAfrica)
    return demog


def find_missing():
    """(xTLH, seed) pairs of the planned 20 runs that have no saved state in the first experiment."""
    import glob, json
    planned = {(x, s) for x in XTLH_VALUES for s in range(NUM_SEEDS)}
    done = set()
    for cfg in glob.glob(os.path.join(manifest.job_directory, OLD_EXP_PREFIX + "*", "*", "config.json")):
        simdir = os.path.dirname(cfg)
        if os.path.exists(os.path.join(simdir, "output", "state-%05d.dtk" % (365 * SERIALIZE_YEARS))):
            p = json.load(open(cfg))["parameters"]
            done.add((round(float(p["x_Temporary_Larval_Habitat"]), 4), int(p["Run_Number"])))
    missing = sorted(planned - done)
    print("Planned runs:", len(planned), " with saved state:", len(done), " missing:", len(missing))
    for m in missing:
        print("  missing: x_Temporary_Larval_Habitat =", m[0], " seed =", m[1])
    return missing


def set_pair(simulation, pair):
    """Sets habitat multiplier and seed together; tags match the first burn-in so the pickup can read both."""
    x, s = pair
    simulation.task.set_parameter("x_Temporary_Larval_Habitat", x)
    simulation.task.set_parameter("Run_Number", s)
    return {"x_Temporary_Larval_Habitat": x, "Run_Number": s}


def general_sim(selected_platform):
    platform = Platform(selected_platform, job_directory=manifest.job_directory,
                        partition=manifest.partition, time=RERUN_TIME,
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
        plugin_report=None
    )
    task.set_sif(manifest.SIF_PATH, platform)

    # climate folder (contains LAKE/2010001-2019365/...) goes in as an asset
    task.common_assets.add_directory(os.path.join(manifest.input_dir, "climate"), relative_path='climate')

    end_day = sim_years*365
    node_ids = pd.read_csv(NODES_FILE)['node_id'].tolist()

    # births, one counter per node
    for n in node_ids:
        add_report_event_counter(task, manifest, start_day=1, end_day=end_day,
                                 node_ids=[n], min_age_years=0, max_age_years=100,
                                 event_trigger_list=["Births"],
                                 filename_suffix="_".join(("node", str(n))))

    # monthly summary (all nodes pooled)
    add_malaria_summary_report(task, manifest, start_day=1, end_day=end_day,
                               reporting_interval=30, age_bins=[0.25, 5, 115],
                               max_number_reports=SERIALIZE_YEARS * 13, filename_suffix='monthly',
                               pretty_format=True)

    # daily spatial report, every node (needed for county-level calibration)
    add_spatial_report_malaria_filtered(task, manifest, start_day=1, end_day=end_day,
                                        reporting_interval=1, node_ids=None,
                                        min_age_years=0.25, max_age_years=100,
                                        spatial_output_channels=["Population",
                                                                 "PCR_Parasite_Prevalence",
                                                                 "New_Clinical_Cases"],
                                        filename_suffix="all_ages")

    user = os.getlogin()
    builder = SimulationBuilder()
    missing = find_missing()
    if not missing:
        raise SystemExit("Nothing missing: all 20 runs already have a saved state. Go to run_lake_pickup3.py.")
    builder.add_sweep_definition(set_pair, missing)
    experiment = Experiment.from_builder(builder, task, name=f'{user}_LAKE_burnin_cal_ext')
    experiment.run(wait_until_done=False, platform=platform)
    print("\n=== RERUN SUBMITTED (missing runs only) ===")
    print("Rerun burn-in ID:", experiment.uid, " simulations:", len(missing))
    print("Write this ID down: run_lake_pickup3.py needs it (BURNIN_RERUN_ID), together with 10b92b9f...")


if __name__ == "__main__":
    import emod_malaria.bootstrap as dtk
    ensure_inputs()
    dtk.setup(pathlib.Path(manifest.eradication_path).parent)
    os.chmod(manifest.eradication_path, 0o755)
    general_sim("SLURM_LOCAL")