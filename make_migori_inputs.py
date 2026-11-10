"""
make_migori_inputs.py - cut the Migori node (node_id 3) out of the 8-node climate files and write a one-node demographics file.
Run once from FE-2026-examples:   python make_migori_inputs.py
Reads  : inputs/climate/LAKE/2010001-2019365/*.bin and *.bin.json   (not changed)
Writes : inputs/climate/LAKE_MIGORI/2010001-2019365/*.bin and *.bin.json
         inputs/demographics/lake_nodes_migori.csv
"""
import os, json
import numpy as np, pandas as pd

NODE, LAT, LON, POP = 3, -0.9995, 34.4744, 5000
SRC = os.path.join("inputs", "climate", "LAKE", "2010001-2019365")
DST = os.path.join("inputs", "climate", "LAKE_MIGORI", "2010001-2019365")
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

dem = os.path.join("inputs", "demographics", "lake_nodes_migori.csv")
pd.DataFrame({"node_id": [NODE], "lat": [LAT], "lon": [LON], "pop": [POP]}).to_csv(dem, index=False)
print("Wrote", dem)
