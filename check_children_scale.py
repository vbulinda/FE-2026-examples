"""
check_children_scale.py - is the spatial 'children' prevalence divided by the children or by everybody?
Run from FE-2026-examples:   python check_children_scale.py
Compares, for ONE pickup simulation and its first 12 months:
  (a) summary report PfPR in the 0.25-5 year bin (properly normalised)
  (b) spatial 'children' prevalence as reported, and divided by the children's share of the population
"""
import glob, json, numpy as np
d = sorted(glob.glob("experiments/e_bulinda_LAKE_pickup_cal_08dab06f*/*/output"))[0]
cfg = json.load(open(d.replace("/output", "/config.json")))["parameters"]
print("simulation", d.split("/")[-2], "xTLH", cfg["x_Temporary_Larval_Habitat"], "seed", cfg["Run_Number"])
j = json.load(open(f"{d}/MalariaSummaryReport_monthly.json"))
b = j["DataByTimeAndAgeBins"]
pf = np.array(b["PfPR by Age Bin"]); pop = np.array(b["Average Population by Age Bin"])
print("age bins:", j["Metadata"]["Age Bins"])
share = (pop[:, 1] - pop[:, 0]) / pop.sum(axis=1)       # 0.5-5 years: bin 0.25-5 minus an estimate of 0.25-0.5 (about the <0.25 bin)
print("children 0.5-5 y share of the population:", share[:12].mean().round(4), "  (0.25-5 y bin share:", (pop[:12, 1] / pop[:12].sum(axis=1)).mean().round(4), ")")
print("(a) summary PfPR, 0.25-5 y, first 12 months, mean: ", pf[:12, 1].mean().round(4))
def rd(n, c):
    a = np.fromfile(f"{d}/SpatialReportMalariaFiltered_{n}_{c}.bin", dtype="<f4"); return a[-4014 * 8:].reshape(4014, 8)
ch = rd("children", "PCR_Parasite_Prevalence")[:365].mean(); al = rd("all_ages", "PCR_Parasite_Prevalence")[:365].mean()
print("(b) spatial children prevalence as reported:        ", ch.round(4), "  all ages:", al.round(4))
print("    children divided by their population share:      ", (ch / share[:12].mean()).round(4))
print("If (a) is close to the last number, the spatial children value is divided by everybody and must be divided by the share.")
print("(a) is a different detection method from PCR, so it is expected to be somewhat lower than the PCR-based number.")
