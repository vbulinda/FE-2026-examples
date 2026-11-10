import os, pandas as pd
d = os.path.join("inputs", "interventions")
t = pd.read_csv(os.path.join(d, "case_management.csv"))
a10, a15 = 0.24, 0.515
for y in range(2010, 2015):
    u5 = a10 + (a15 - a10) * (y - 2010) / 5
    for lo, f in [(0, 1.0), (5, 0.6), (15, 0.4)]:
        m = (t.year == y) & (t.trigger == "NewClinicalCase") & (t.age_min == lo)
        t.loc[m, "coverage"] = round(u5 * f, 4)
        t.loc[m, "source"] = "u5 ACT: Kenya MIS 2010 lake endemic 24% (2010) rising linearly to MIS 2015 value (assumption); older ages x0.6 / x0.4 (assumption)"
t.to_csv(os.path.join(d, "case_management_v2.csv"), index=False)
print(t[(t.node_id == 1) & (t.trigger == "NewClinicalCase")].groupby(["year", "age_min"]).coverage.first().unstack().head(7))
