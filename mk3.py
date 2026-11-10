s = open("run_bungoma_pickup_v2.py").read()
old = '''    for y, sd in ((2021, 4180), (2024, 5291)):
        n_df = pd.concat([n_df, base.assign(year=y, start_day=sd, coverage=0.825)],
                         ignore_index=True)
'''
new = '''    NET_EXTRA = [int(v) for v in os.environ.get("NET_EXTRA", "").split(",") if v]
    n_df.loc[n_df["year"] == 2017, "start_day"] = int(os.environ.get("NET_2017", "2800"))
    extra = [(2010 + v // 365, v, 0.592) for v in NET_EXTRA]
    for y, sd, cv in [(2021, 4180, 0.825), (2024, 5291, 0.825)] + extra:
        n_df = pd.concat([n_df, base.assign(year=y, start_day=sd, coverage=cv)],
                         ignore_index=True)
    print("ITN rounds (year, start_day, coverage):")
    print(n_df[["year", "start_day", "coverage"]].to_string())
'''
assert s.count(old) == 1
s = s.replace(old, new)
a = 'if y <= year] or [year]'
assert s.count(a) == 1
s = s.replace(a, 'if y <= year] or [int(itn_age["year"].min())]')
b = '_pickup_v2_nasa_{SCENARIO}_imp{IMPORT_PER_MONTH}'
assert s.count(b) == 1
s = s.replace(b, '_pickup_v3_nasa_{SCENARIO}_imp{IMPORT_PER_MONTH}_net{NET_TAG}')
s = s.replace('SCENARIO = os.environ.get(', 'NET_TAG = os.environ.get("NET_TAG", "base")\nSCENARIO = os.environ.get(', 1)
open("run_bungoma_pickup_v3.py", "w").write(s)
print("wrote run_bungoma_pickup_v3.py", len(s.splitlines()), "lines")
