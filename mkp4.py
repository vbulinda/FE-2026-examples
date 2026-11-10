# mkp4.py - run_migori_pickup_v3.py -> run_migori_pickup_v4.py : env IRS_KILL overrides killing_initial_effect (table value 0.8) for every IRS round.
s = open("run_migori_pickup_v3.py").read()
def rep(a, b):
    global s
    assert s.count(a) == 1, a
    s = s.replace(a, b)
rep('IRS_DAYS = int(os.environ.get("IRS_DAYS", "0"))',
    'IRS_KILL = float(os.environ.get("IRS_KILL", "0"))   # 0 = table kill_effect (0.8); else this killing_initial_effect for all rounds\n'
    'IRS_DAYS = int(os.environ.get("IRS_DAYS", "0"))')
rep('killing_initial_effect=float(n_df["kill_effect"][r]),', 'killing_initial_effect=(IRS_KILL or float(n_df["kill_effect"][r])),')
rep('full-effect days {dur}")', 'full-effect days {dur}, kill effect {IRS_KILL or float(n_df[\'kill_effect\'][r])}")')
open("run_migori_pickup_v4.py", "w").write(s)
print("wrote run_migori_pickup_v4.py")
