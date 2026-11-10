"""
make_bungoma_inputs_nasa.py - cut the siaya node (node_id 2) out of the 8-node climate files and write a one-node demographics file.
Run once from FE-2026-examples:   python make_bungoma_inputs_nasa.py
Reads  : inputs/climate/LAKE_NASA/2010001-2026304/*.bin and *.bin.json   (not changed)
Writes : inputs/climate/LAKE_SIAYA_NASA/2010001-2026304/*.bin and *.bin.json
         inputs/demographics/lake_nodes_siaya.csv
"""
import os, json
import numpy as np, pandas as pd

NODE, LAT, LON, POP = 2, 0.0611, 34.2882, 5000      # siaya, from COUNTY_NODES; POP = TEST_POP of the calibration runs
SRC = os.path.join("inputs", "climate", "LAKE_NASA", "2010001-2026304")
DST = os.path.join("inputs", "climate", "LAKE_SIAYA_NASA", "2010001-2026304")
STEMS = ["dtk_15arcmin_air_temperature_daily", "dtk_15arcmin_land_temperature_daily",
         "dtk_15arcmin_relative_humidity_daily", "dtk_15arcmin_rainfall_daily"]
os.makedirs(DST, exist_ok=True)

for stem in STEMS:
    meta = json.load(open(os.path.join(SRC, stem + ".bin.json")))
    m = meta["Metadata"]; days, n = m["DatavalueCount"], m["NodeCount"]
    ids = [int(meta["NodeOffsets"][16 * i:16 * i + 8], 16) for i in range(n)]
    offs = [int(meta["NodeOffsets"][16 * i + 8:16 * i + 16], 16) for i in range(n)]
    raw = np.fromfile(os.path.join(SRC, stem + ".bin"), dtype="<f4")
    assert raw.size == n * days, f"{stem}: expected {n*days} values, found {raw.size}"
    k = ids.index(NODE)
    arr = raw[offs[k] // 4: offs[k] // 4 + days].astype("<f4")
    arr.tofile(os.path.join(DST, stem + ".bin"))
    m2 = dict(m)
    m2.update({"WeatherCellCount": 1, "NodeCount": 1, "OffsetEntryCount": 1, "NumberDTKNodes": 1,
               "UpperLatitude": LAT, "BottomLatitude": LAT, "LeftLongitude": LON, "RightLongitude": LON})
    with open(os.path.join(DST, stem + ".bin.json"), "w") as f:
        json.dump({"Metadata": m2, "NodeOffsets": f"{NODE:08x}{0:08x}"}, f, indent=4)
    print(f"{stem}: node {NODE} -> {arr.size} days, range {arr.min():.2f} to {arr.max():.2f}")

dem = os.path.join("inputs", "demographics", "lake_nodes_siaya.csv")
pd.DataFrame({"node_id": [NODE], "lat": [LAT], "lon": [LON], "pop": [POP]}).to_csv(dem, index=False)
print("Wrote", dem)
