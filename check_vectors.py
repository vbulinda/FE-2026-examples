import glob, json
import pandas as pd
for exp in glob.glob("experiments/*pickup_v1_036ac7c8*"):
    for sim in glob.glob(exp + "/*/config.json"):
        c = json.load(open(sim))["parameters"]
        if abs(c["x_Temporary_Larval_Habitat"] - 8.9797) < 0.01 and c["Run_Number"] == 0:
            out = sim.replace("config.json", "output/")
            v = pd.read_csv(out + "ReportVectorStats.csv")
            v["year"] = 2010 + (v["Time"] - 1) // 365
            print("xTLH 8.9797, seed 0 | yearly mean of daily values")
            print(v.groupby(["year", "Species"])[["VectorPopulation", "STATE_INFECTIOUS"]].mean().round(1).unstack().to_string())
            m = json.load(open(out + "MalariaSummaryReport_monthly.json"))["DataByTime"]
            e = pd.DataFrame({"t": m["Time Of Report"], "eir": m["Annual EIR"]})
            e["year"] = 2010 + (e["t"] - 1) // 365
            print("\nAnnual EIR, yearly mean of monthly reports")
            print(e.groupby("year")["eir"].mean().round(2).to_string())
