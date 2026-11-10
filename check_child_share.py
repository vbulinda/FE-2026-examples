import glob, json
import numpy as np
sim = glob.glob("experiments/*pickup_v1_036ac7c8*/*/output/MalariaSummaryReport_monthly.json")[0]
j = json.load(open(sim))
print("top-level keys:", list(j.keys()))
d = j["DataByTimeAndAgeBins"]
print("age-bin keys:", list(d.keys()))
pop = np.array(d["Average Population by Age Bin"], dtype=float)
print("age-bin array shape:", pop.shape)
print("share aged 0.25-5 y of the node (mean over months): %.4f" % np.mean(pop[:, 0] / pop.sum(axis=1)))
dem = glob.glob("experiments/*pickup_v1_036ac7c8*/**/demographics.json", recursive=True)
if dem:
    a = json.load(open(dem[0]))["Defaults"]["IndividualAttributes"]["AgeDistribution"]
    ages = np.array(a["ResultValues"], dtype=float)
    cum = np.array(a["DistributionValues"], dtype=float)
    print("age distribution file:", dem[0])
    print("share aged 0.5-5 y at start: %.4f" % (np.interp(5 * 365, ages, cum) - np.interp(0.5 * 365, ages, cum)))
else:
    print("no demographics.json found")
