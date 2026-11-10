s = open("run_migori_structure3.py").read()


def rep(old, new, n=1):
    global s
    assert s.count(old) == n, ("NOT FOUND or wrong count:", old[:70], s.count(old))
    s = s.replace(old, new)


i0 = s.index('"""')
i1 = s.index('"""', i0 + 3) + 3
DOC = ('"""\nrun_bungoma_pickup_v2.py - Bungoma (node 7) pickup, 6,144 days (2010 to\n'
       '31 Oct 2026), real NASA POWER weather, no IRS. Built from\n'
       'run_migori_structure3.py; the changes are in the calibration log,\n'
       'steps 45 and 51.   Run: python run_bungoma_pickup_v2.py\n"""')
s = DOC + s[i1:]
rep('NODE = 3\n', 'NODE = 7\n')
rep('BURNIN_EXP_IDS = ["d44be5ca-8222-43dc-a98d-822b23ae8157"]',
    'BURNIN_EXP_IDS = ["b4f301d5-43e7-4ed8-98be-d3b6fac618a2"]')
rep('ARA_MULT = 0.125\n', 'G_MULT, ARA_MULT, FUN_MULT = 20, 0.0225, 0.006   # as in the burn-in\n')
rep('IMPORT_PER_MONTH = int(os.environ.get("IMPORT", "0"))',
    'IMPORT_PER_MONTH = int(os.environ.get("IMPORT", "1"))')
rep('LAKE_MIGORI", "2010001-2019365"', 'LAKE_BUNGOMA_NASA", "2010001-2026304"')
rep('lake_nodes_migori.csv', 'lake_nodes_bungoma.csv')
rep('lake_nodes_migori_test.csv', 'lake_nodes_bungoma_test.csv', 2)
SP_OLD = ('    conf.set_max_larval_capacity(config, "arabiensis", "TEMPORARY_RAINFALL",'
          ' int(800000000 * ARA_MULT))\n'
          '    conf.set_max_larval_capacity(config, "arabiensis", "CONSTANT",'
          ' int(80000000 * ARA_MULT))\n')
SP_NEW = ('    sp = conf.set_max_larval_capacity\n'
          '    sp(config, "gambiae", "WATER_VEGETATION", int(20000000 * G_MULT))\n'
          '    sp(config, "arabiensis", "TEMPORARY_RAINFALL", int(800000000 * ARA_MULT))\n'
          '    sp(config, "arabiensis", "CONSTANT", int(80000000 * ARA_MULT))\n'
          '    sp(config, "funestus", "TEMPORARY_RAINFALL", int(800000000 * FUN_MULT))\n'
          '    sp(config, "funestus", "CONSTANT", int(80000000 * FUN_MULT))\n')
rep(SP_OLD, SP_NEW)
rep('    if os.path.exists(irs_path) and SCENARIO not in ("no_irs", "no_irs_itn17"):',
    '    if False:   # Bungoma has no IRS (no node 7 rows in IRS.json)')
rep('print("NOTE: no IRS.json found: this run has NO IRS.")',
    'print("NOTE: Bungoma has no IRS rows (node 7): this run has NO IRS.")')
rep('Keep only the Migori rows.', 'Keep only the rows of NODE.')
rep('PICKUP_YEARS = 11', 'SIM_DAYS = 6144   # 1 Jan 2010 to 31 Oct 2026; equals the climate days')
rep('config.parameters.Simulation_Duration = PICKUP_YEARS * 365 + 30',
    'config.parameters.Simulation_Duration = SIM_DAYS')
rep('repetitions=PICKUP_YEARS * 365 // 30', 'repetitions=SIM_DAYS // 30')
rep('end_day = PICKUP_YEARS * 365', 'end_day = SIM_DAYS')
rep('max_number_reports=PICKUP_YEARS * 13 + 1', 'max_number_reports=SIM_DAYS // 30 + 2')
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
rep('_LAKE_MIGORI_pickup_v1_{SCENARIO}_irs{int(IRS_REAL)}_imp{IMPORT_PER_MONTH}'
    '_d{IRS_DAYS}_c{IRS_2017_COVERAGE}',
    '_LAKE_BUNGOMA_pickup_v2_nasa_{SCENARIO}_imp{IMPORT_PER_MONTH}')
rep('=== MIGORI PICKUP SUBMITTED ===', '=== BUNGOMA PICKUP SUBMITTED ===')
open("run_bungoma_pickup_v2.py", "w").write(s)
print("written")
