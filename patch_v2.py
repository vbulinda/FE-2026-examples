s = open("run_bungoma_pickup.py").read()


def rep(old, new, n=1):
    global s
    assert s.count(old) == n, ("NOT FOUND or wrong count:", old[:70], s.count(old))
    s = s.replace(old, new)


i0 = s.index('"""')
i1 = s.index('"""', i0 + 3) + 3
DOC = ('"""\nrun_bungoma_pickup_v2.py - Bungoma pickup, 6,144 days (2010 to 31 Oct 2026),\n'
       'real NASA POWER weather. Built from run_bungoma_pickup.py; the changes\n'
       'are listed in the calibration log, step 51.\n'
       'Run: python run_bungoma_pickup_v2.py\n"""')
s = DOC + s[i1:]
rep('PICKUP_YEARS = 11', 'SIM_DAYS = 6144   # 1 Jan 2010 to 31 Oct 2026; equals the climate days')
rep('config.parameters.Simulation_Duration = PICKUP_YEARS * 365 + 30',
    'config.parameters.Simulation_Duration = SIM_DAYS')
rep('repetitions=PICKUP_YEARS * 365 // 30', 'repetitions=SIM_DAYS // 30')
rep('end_day = PICKUP_YEARS * 365', 'end_day = SIM_DAYS')
rep('max_number_reports=PICKUP_YEARS * 13 + 1', 'max_number_reports=SIM_DAYS // 30 + 2')
rep('LAKE_BUNGOMA", "2010001-2019365"', 'LAKE_BUNGOMA_NASA", "2010001-2026304"')
CM_OLD = '    n_df = cm_all[cm_all["phase"] == "pickup"].reset_index(drop=True)\n'
CM_NEW = (CM_OLD +
          '    ext = []\n'
          '    for y in range(2021, 2027):   # 2020 rows carried forward, log step 51\n'
          '        rows = n_df[n_df["year"] == 2020].copy()\n'
          '        rows["year"] = y\n'
          '        rows["start_day"] = (y - 2010) * 365 + 1\n'
          '        ext.append(rows)\n'
          '    n_df = pd.concat([n_df] + ext, ignore_index=True)\n')
rep(CM_OLD, CM_NEW)
IT_OLD = '    n_df = itn_all[itn_all["phase"] == "pickup"].reset_index(drop=True)\n'
IT_NEW = (IT_OLD +
          '    base = n_df[n_df["year"] == 2017].iloc[[0]]   # PBO benefit not modelled\n'
          '    for y, sd in ((2021, 4180), (2024, 5291)):\n'
          '        n_df = pd.concat([n_df, base.assign(year=y, start_day=sd, coverage=0.825)],\n'
          '                         ignore_index=True)\n')
rep(IT_OLD, IT_NEW)
rep('age_y = itn_age[itn_age["year"] == year]',
    'yrs = [y for y in itn_age["year"].unique() if y <= year] or [year]\n'
    '        age_y = itn_age[itn_age["year"] == max(yrs)]')
rep('_LAKE_BUNGOMA_pickup_v1_{SCENARIO}_imp{IMPORT_PER_MONTH}',
    '_LAKE_BUNGOMA_pickup_v2_nasa_{SCENARIO}_imp{IMPORT_PER_MONTH}')
open("run_bungoma_pickup_v2.py", "w").write(s)
print("written")
