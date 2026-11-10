import glob, json
import numpy as np
sims = sorted(glob.glob("experiments/*pickup_v1_036ac7c8*/*/output/MalariaSummaryReport_monthly.json"))
shares = []
for s in sims:
    pop = np.array(json.load(open(s))["DataByTimeAndAgeBins"]["Average Population by Age Bin"], dtype=float)
    shares.append(pop[:, 1] / pop.sum(axis=1))
sh = np.array(shares)
print("simulations:", len(sims), "| months:", sh.shape[1])
print("model share aged 0.25-5 y, mean of all months and runs: %.4f" % sh.mean())
print("model share aged 0.25-5 y, months 66-68 (the Jul-Aug 2015 window): %.4f" % sh[:, 66:69].mean())
print("model share aged 0.25-5 y, last 3 months (2020): %.4f" % sh[:, -4:-1].mean())
dem = glob.glob("experiments/*pickup_v1_036ac7c8*/**/demographics.json", recursive=True)[0]
a = json.load(open(dem))["Defaults"]["IndividualAttributes"]["AgeDistribution"]
ages = np.array(a["ResultValues"], dtype=float); cum = np.array(a["DistributionValues"], dtype=float)
sh_025 = np.interp(5 * 365, ages, cum) - np.interp(0.25 * 365, ages, cum)
sh_05 = np.interp(5 * 365, ages, cum) - np.interp(0.5 * 365, ages, cum)
print("initial distribution: 0.25-5 y %.4f | 0.5-5 y %.4f | ratio %.4f" % (sh_025, sh_05, sh_05 / sh_025))
print("model share aged 0.5-5 y (estimate = model 0.25-5 x ratio): %.4f" % (sh[:, 66:69].mean() * sh_05 / sh_025))
