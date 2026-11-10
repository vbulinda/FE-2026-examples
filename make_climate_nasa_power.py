"""
make_climate_nasa_power.py - rebuild the daily climate files for all 8 nodes with real NASA POWER weather,
2010-01-01 to 2026-10-31 (6,144 days: 365-day years, 29 February removed), so that runs no longer repeat 2010 weather in 2020.
Same recipe as the existing files (inputs/climate/LAKE/2010001-2019365: NASA POWER daily point data T2M, RH2M, PRECTOTCORR):
  air and land temperature = T2M (C); relative humidity = RH2M / 100 (fraction); rainfall = PRECTOTCORR (mm/day).
Step 1: downloads each node (lat, lon from inputs/demographics/lake_nodes.csv) and keeps the raw JSON in inputs/climate/LAKE_NASA/raw_nasa_power/.
Step 2: checks 2010-2019 of the new series against the existing files (two leap-day rules are tried; the better one is used).
Step 3: writes the new files to inputs/climate/LAKE_NASA/2010001-2026304/ (nothing existing is changed).
Run once from FE-2026-examples (emod-venv active):   python make_climate_nasa_power.py
"""
import os, json, time, datetime, urllib.request
import numpy as np, pandas as pd

START, END, N_DAYS = "20100101", "20261031", 6144
NODES_CSV = os.path.join("inputs", "demographics", "lake_nodes.csv")
OLD = os.path.join("inputs", "climate", "LAKE", "2010001-2019365")
BASE = os.path.join("inputs", "climate", "LAKE_NASA")
NEW, RAW = os.path.join(BASE, "2010001-2026304"), os.path.join(BASE, "raw_nasa_power")
STEMS = {"air": "dtk_15arcmin_air_temperature_daily", "land": "dtk_15arcmin_land_temperature_daily",
         "rh": "dtk_15arcmin_relative_humidity_daily", "rain": "dtk_15arcmin_rainfall_daily"}
URL = ("https://power.larc.nasa.gov/api/temporal/daily/point?parameters=T2M,RH2M,PRECTOTCORR&community=AG"
       "&longitude={lon}&latitude={lat}&start=" + START + "&end=" + END + "&format=JSON")


def download(node, lat, lon):
    path = os.path.join(RAW, f"node_{node}.json")
    if not os.path.exists(path):
        for attempt in range(4):
            try:
                with urllib.request.urlopen(URL.format(lat=lat, lon=lon), timeout=180) as r:
                    open(path, "wb").write(r.read())
                break
            except Exception as e:
                print("  retry", attempt + 1, e); time.sleep(5)
        time.sleep(1)
    return json.load(open(path))["properties"]["parameter"]


def to_arrays(par, rule):
    """rule 'feb29' drops 29 February; rule 'dec31' drops 31 December of leap years. Returns dict of 6,144-value arrays."""
    dates = sorted(par["T2M"])
    leap_dec31 = {f"{y}1231" for y in range(2010, 2027) if y % 4 == 0}
    keep = [d for d in dates if (d[4:] != "0229" if rule == "feb29" else d not in leap_dec31)]
    out, nfill = {}, 0
    for k in ("T2M", "RH2M", "PRECTOTCORR"):
        s = pd.Series([par[k][d] for d in keep], dtype=float)
        bad = (s <= -900); nfill += int(bad.sum())
        out[k] = s.mask(bad).interpolate(limit_direction="both").values
    assert len(keep) == N_DAYS, (rule, len(keep))
    return {"air": out["T2M"], "land": out["T2M"], "rh": np.clip(out["RH2M"] / 100.0, 0, 1), "rain": out["PRECTOTCORR"]}, nfill


def read_old(stem):
    meta = json.load(open(os.path.join(OLD, stem + ".bin.json")))
    days, n = meta["Metadata"]["DatavalueCount"], meta["Metadata"]["NodeCount"]
    ids = [int(meta["NodeOffsets"][16 * i:16 * i + 8], 16) for i in range(n)]
    offs = [int(meta["NodeOffsets"][16 * i + 8:16 * i + 16], 16) // 4 for i in range(n)]
    raw = np.fromfile(os.path.join(OLD, stem + ".bin"), dtype="<f4")
    return meta, {i: raw[o:o + days].astype(float) for i, o in zip(ids, offs)}, ids


if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True); os.makedirs(NEW, exist_ok=True)
    nodes = pd.read_csv(NODES_CSV)
    old = {k: read_old(v) for k, v in STEMS.items()}
    ids = old["air"][2]
    series = {}
    for _, r in nodes.iterrows():
        print(f"node {int(r.node_id)}: lat {r.lat}, lon {r.lon}")
        series[int(r.node_id)] = download(int(r.node_id), r.lat, r.lon)
    # ---- validation against the existing 2010-2019 files, two leap-day rules ----
    scores = {}
    for rule in ("feb29", "dec31"):
        tot = 0.0
        print(f"\nLeap-day rule '{rule}': mean absolute difference, new minus existing, 2010-2019 (days 1-3650)")
        print("node | air C | land C | RH frac | rain mm/day | rain correlation")
        for nid in ids:
            arr, nfill = to_arrays(series[nid], rule)
            row = []
            for k in ("air", "land", "rh", "rain"):
                a, b = arr[k][:3650], old[k][1][nid]
                row.append(float(np.mean(np.abs(a - b)))); tot += row[-1] / max(float(np.std(b)), 1e-6)
            cor = float(np.corrcoef(arr["rain"][:3650], old["rain"][1][nid])[0, 1])
            print(f"{nid:4d} | {row[0]:.3f} | {row[1]:.3f} | {row[2]:.4f} | {row[3]:.3f} | {cor:.3f}   (filled values: {nfill})")
        scores[rule] = tot
    rule = min(scores, key=scores.get)
    print("\nLeap-day rule used:", rule, "| scores (lower = closer to the existing files):", {k: round(v, 2) for k, v in scores.items()})
    # ---- write the new files ----
    now = datetime.date.today().isoformat()
    for k, stem in STEMS.items():
        meta, _, _ = old[k]
        m = dict(meta["Metadata"])
        m.update({"DatavalueCount": N_DAYS, "OriginalDataYears": "2010-2026 (to 31 October 2026)",
                  "DataProvenance": f"NASA POWER daily point data (T2M, RH2M, PRECTOTCORR), API v2.10.3, downloaded {now}; leap-day rule {rule}",
                  "DateCreated": now, "Author": "Lake Victoria EMOD project"})
        data = np.concatenate([to_arrays(series[nid], rule)[0][k] for nid in ids]).astype("<f4")
        data.tofile(os.path.join(NEW, stem + ".bin"))
        offsets = "".join(f"{nid:08x}{i * N_DAYS * 4:08x}" for i, nid in enumerate(ids))
        json.dump({"Metadata": m, "NodeOffsets": offsets}, open(os.path.join(NEW, stem + ".bin.json"), "w"), indent=4)
        print(f"{stem}: {len(ids)} nodes x {N_DAYS} days, range {data.min():.2f} to {data.max():.2f}")
    print("\nWritten to", NEW)
