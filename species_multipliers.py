"""
species_multipliers.py - habitat multipliers for a target species mix, from the two vector-only runs of run_<county>_species_unit.py.
Run from FE-2026-examples after both jobs finish:   python species_multipliers.py kisumu
Finds the two experiment folders experiments/e_*LAKE_<COUNTY>_species_* (names are truncated by idmtools, so the runs are told apart by their arabiensis density) and reads each ReportVectorStats.csv (Time >= 730; vectors per day = summed
VectorPopulation / number of days, per species). Vectors per day scale in direct proportion to each species' multiplier
(Homa Bay species grid, log Steps 121-122), so from the unit run (g 1, ara 1, fun 1) the multipliers for any share s_i are
m_i = s_i x T / u_i, where T is the total vectors per day of the Migori multipliers (g 1, ara 0.125, fun 1) at THIS county's
climate (as in Step 122, so xTLH stays comparable). The second run checks linearity: predicted against observed total.
These are model scenarios, not field measurements. Prints candidate mixes; the choice is made afterwards.
"""
import glob, sys
import pandas as pd

county = sys.argv[1].lower()
SP = ["gambiae", "arabiensis", "funestus"]
U = county.upper()
dirs = sorted(glob.glob(f"experiments/e_*LAKE_{U}_species_*"))
if len(dirs) != 2:
    raise SystemExit(f"expected exactly two experiment folders experiments/e_*LAKE_{U}_species_*, found {len(dirs)}: {dirs}")
dens = []
for d in dirs:
    fs = glob.glob(d + "/**/ReportVectorStats.csv", recursive=True)
    if not fs:
        raise SystemExit(f"{d}: not finished or ReportVectorStats.csv not found")
    print("read", d)
    df = pd.read_csv(fs[0]); df = df[df["Time"] >= 730]
    days = df["Time"].nunique()
    by = df.groupby("Species")["VectorPopulation"].sum()
    dens.append({s: float(by.get(s, 0.0)) / days for s in SP})
# the two runs are told apart by their data: the Migori-multiplier run (ara 0.125) has 8 times fewer arabiensis than the unit run (ara 1)
dens.sort(key=lambda x: -x["arabiensis"])
r = dens[0]["arabiensis"] / dens[1]["arabiensis"]
if not 7.2 <= r <= 8.8:
    raise SystemExit(f"arabiensis ratio between the two runs is {r:.2f}, expected about 8: these are not the unit and Migori-multiplier runs")
runs = {(1.0, 1.0, 1.0): dens[0], (1.0, 0.125, 1.0): dens[1]}
unit, migori = runs[(1.0, 1.0, 1.0)], runs[(1.0, 0.125, 1.0)]
print(f"{county}: vectors per day per unit multiplier (gambiae, arabiensis, funestus):",
      {s: round(unit[s]) for s in SP})
pred = unit["gambiae"] * 1 + unit["arabiensis"] * 0.125 + unit["funestus"] * 1
T = sum(migori.values())
print(f"Migori multipliers (g 1, ara 0.125, fun 1): observed total {T:,.0f} per day, predicted from the unit run {pred:,.0f} "
      f"(difference {100 * (T - pred) / pred:+.1f} percent; linearity holds if within a few percent)")
print("Migori multipliers give shares (gam, ara, fun) %:", tuple(round(100 * migori[s] / T, 1) for s in SP))
MIXES = {"Kisian 2018 (Kisumu): gam 32, ara 49, fun 19": (0.32, 0.49, 0.19),
         "Bungoma-like: gam 97.2, ara 2.2, fun 0.6": (0.972, 0.022, 0.006),
         "Funestus-dominant (Migori-like): gam 5.3, ara 10.5, fun 84.2": (0.053, 0.105, 0.842)}
print("\nMultipliers that give each mix at total density T (use the rounded values in the burn-in):")
for name, sh in MIXES.items():
    m = [sh[i] * T / unit[SP[i]] for i in range(3)]
    print(f"  {name}:  G_MULT {m[0]:.4g}  ARA_MULT {m[1]:.4g}  FUN_MULT {m[2]:.4g}")
