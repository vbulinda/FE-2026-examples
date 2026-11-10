# mkp3.py - run_migori_pickup_v2.py -> run_migori_pickup_v3.py : per-round IRS coverage and full-effect days.
# Env IRS_COVS="0.977,0.83,0.83,0.904" and IRS_DAYS_LIST="300,300,300,150" (one value per round, in date order; need IRS_REAL=1 = 4 rounds).
s = open("run_migori_pickup_v2.py").read()
def rep(a, b):
    global s
    assert s.count(a) == 1, a
    s = s.replace(a, b)
rep('IRS_DAYS = int(os.environ.get("IRS_DAYS", "0"))',
    'IRS_COVS = [float(x) for x in os.environ.get("IRS_COVS", "").split(",") if x]   # per-round coverage (overrides IRS_COV)\n'
    'IRS_DAYS_LIST = [int(x) for x in os.environ.get("IRS_DAYS_LIST", "").split(",") if x]   # per-round full-effect days (overrides IRS_DAYS)\n'
    'IRS_DAYS = int(os.environ.get("IRS_DAYS", "0"))')
rep('            print(f"IRS node {NODE}: start_day {sd}, coverage used {cov} (table value {float(n_df[\'coverage\'][r])}), full-effect days {int(IRS_DAYS or n_df[\'kill_duration\'][r])}")',
    '            if IRS_COVS:\n                assert len(IRS_COVS) == len(n_df), "IRS_COVS needs one value per round"\n                cov = IRS_COVS[r]\n'
    '            dur = IRS_DAYS_LIST[r] if IRS_DAYS_LIST else int(IRS_DAYS or n_df["kill_duration"][r])\n'
    '            if IRS_DAYS_LIST:\n                assert len(IRS_DAYS_LIST) == len(n_df), "IRS_DAYS_LIST needs one value per round"\n'
    '            print(f"IRS node {NODE}: start_day {sd}, coverage used {cov} (table value {float(n_df[\'coverage\'][r])}), full-effect days {dur}")')
rep('killing_box_duration=int(IRS_DAYS or n_df["kill_duration"][r]),', 'killing_box_duration=dur,')
open("run_migori_pickup_v3.py", "w").write(s)
print("wrote run_migori_pickup_v3.py")
