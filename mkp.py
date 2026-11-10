s = open("run_migori_structure3.py").read()


def rep(old, new, n=1):
    global s
    assert s.count(old) == n, ("NOT FOUND or wrong count:", old[:70], s.count(old))
    s = s.replace(old, new)


i1 = s.index('"""', s.index('"""') + 3) + 3
DOC = ('"""\nrun_migori_pickup_v2.py - Migori (node 3) pickup, 6,144 days (2010 to 31 Oct 2026),\n'
       'real NASA POWER weather, burn-in v2 (baseline case management + importation).\n'
       'Env: BURNIN_ID (required), IRS_COV, IRS_DAYS, IRS_REAL, IMPORT, NET_EXTRA, NET_2017,\n'
       'NET_LATE, ITN_KAPPA, ITN_LAMBDA, NET_TAG.   Run: BURNIN_ID=<id> python run_migori_pickup_v2.py\n"""')
s = DOC + s[i1:]
rep('BURNIN_EXP_IDS = ["d44be5ca-8222-43dc-a98d-822b23ae8157"]',
    'BURNIN_EXP_IDS = [os.environ.get("BURNIN_ID", "")]\n'
    'if not BURNIN_EXP_IDS[0]:\n    raise SystemExit("Set BURNIN_ID=<burn-in experiment id>")\n'
    'NET_TAG = os.environ.get("NET_TAG", "base")')
rep('IMPORT_PER_MONTH = int(os.environ.get("IMPORT", "0"))',
    'IMPORT_PER_MONTH = int(os.environ.get("IMPORT", "1"))')
rep('LAKE_MIGORI", "2010001-2019365"', 'LAKE_MIGORI_NASA", "2010001-2026304"')
rep('PICKUP_YEARS = 11', 'SIM_DAYS = 6144   # 1 Jan 2010 to 31 Oct 2026; equals the climate days')
rep('config.parameters.Simulation_Duration = PICKUP_YEARS * 365 + 30',
    'config.parameters.Simulation_Duration = SIM_DAYS')
rep('repetitions=PICKUP_YEARS * 365 // 30', 'repetitions=SIM_DAYS // 30')
rep('end_day = PICKUP_YEARS * 365', 'end_day = SIM_DAYS')
rep('max_number_reports=PICKUP_YEARS * 13 + 1', 'max_number_reports=SIM_DAYS // 30 + 2')
CM_OLD = '    n_df = cm_all[cm_all["phase"] == "pickup"].reset_index(drop=True)\n'
CM_NEW = (CM_OLD +
          '    ext = []\n'
          '    for y in range(2021, 2027):   # 2020 rows carried forward\n'
          '        rows = n_df[n_df["year"] == 2020].copy()\n'
          '        rows["year"] = y\n'
          '        rows["start_day"] = (y - 2010) * 365 + 1\n'
          '        ext.append(rows)\n'
          '    n_df = pd.concat([n_df] + ext, ignore_index=True)\n')
rep(CM_OLD, CM_NEW)
IT_OLD = '    n_df = itn_all[itn_all["phase"] == "pickup"].reset_index(drop=True)\n'
IT_NEW = (IT_OLD +
          '    base = n_df[n_df["year"] == 2017].iloc[[0]]\n'
          '    cov0 = float(base["coverage"].iloc[0])\n'
          '    n_df.loc[n_df["year"] == 2017, "start_day"] = int(os.environ.get("NET_2017", "2661"))\n'
          '    NET_EXTRA = [int(v) for v in os.environ.get("NET_EXTRA", "").split(",") if v]\n'
          '    late = [] if os.environ.get("NET_LATE", "1") == "0" else [(2021, 4180), (2024, 5291)]\n'
          '    for y, sd in late + [(2010 + v // 365, v) for v in NET_EXTRA]:\n'
          '        n_df = pd.concat([n_df, base.assign(year=y, start_day=sd, coverage=cov0)],\n'
          '                         ignore_index=True)\n'
          '    print("ITN rounds (year, start_day, coverage):")\n'
          '    print(n_df[["year", "start_day", "coverage"]].to_string())\n')
rep(IT_OLD, IT_NEW)
rep('age_y = itn_age[itn_age["year"] == year]',
    'yrs = [y for y in itn_age["year"].unique() if y <= year] or [int(itn_age["year"].min())]\n'
    '        age_y = itn_age[itn_age["year"] == max(yrs)]')
rep('"Expiration_Period_Kappa": _num(sub["discard_k"][0]),',
    '"Expiration_Period_Kappa": float(os.environ.get("ITN_KAPPA", "2")),')
rep('"Expiration_Period_Lambda": _num(sub["discard_l"][0])}',
    '"Expiration_Period_Lambda": float(os.environ.get("ITN_LAMBDA", "1254"))}')
rep('itn.add_scheduled_usage_dependent_bednet(',
    'print("ITN discard", year, discard)\n        itn.add_scheduled_usage_dependent_bednet(')
rep('_LAKE_MIGORI_pickup_v1_{SCENARIO}', '_LAKE_MIGORI_pickup_v2_{SCENARIO}')
rep('_c{IRS_2017_COVERAGE}', '_c{IRS_2017_COVERAGE}_net{NET_TAG}')
open("run_migori_pickup_v2.py", "w").write(s)
print("wrote run_migori_pickup_v2.py", len(s.splitlines()), "lines")
