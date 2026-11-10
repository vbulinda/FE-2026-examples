import json, os
import numpy as np
D = os.path.join("inputs", "climate", "LAKE_NASA", "2010001-2026304")
def get(stem, node):
    meta = json.load(open(os.path.join(D, stem + ".bin.json")))
    days, n = meta["Metadata"]["DatavalueCount"], meta["Metadata"]["NodeCount"]
    ids = [int(meta["NodeOffsets"][16*i:16*i+8], 16) for i in range(n)]
    offs = [int(meta["NodeOffsets"][16*i+8:16*i+16], 16)//4 for i in range(n)]
    raw = np.fromfile(os.path.join(D, stem + ".bin"), dtype="<f4")
    return raw[offs[ids.index(node)]:offs[ids.index(node)]+days]
for node, name in ((3, "Migori"), (7, "Bungoma")):
    t = get("dtk_15arcmin_air_temperature_daily", node); r = get("dtk_15arcmin_rainfall_daily", node)
    print(f"\n{name} (node {node}): days {len(t)}")
    print("year | mean air C | rain mm/year | days with rain >1 mm")
    for y in range(2010, 2027):
        i0 = (y - 2010) * 365; i1 = min(i0 + 365, len(t))
        print(f"{y} | {t[i0:i1].mean():10.2f} | {r[i0:i1].sum():12.0f} | {(r[i0:i1] > 1).sum():4d}" + ("   (to 31 Oct)" if y == 2026 else ""))
    a, b = r[0:365], r[3650:4015]
    print(f"2010 vs 2020 daily rainfall: correlation {np.corrcoef(a, b)[0,1]:.2f} (1.00 would mean the same weather); max daily rain {r.max():.0f} mm on day {int(r.argmax())} = {2010 + int(r.argmax())//365}")
