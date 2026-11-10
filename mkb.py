s = open("run_migori_burnin.py").read()


def rep(old, new, n=1):
    global s
    assert s.count(old) == n, ("NOT FOUND or wrong count:", old[:70], s.count(old))
    s = s.replace(old, new)


i1 = s.index('"""', s.index('"""') + 3) + 3
DOC = ('"""\nrun_migori_burnin_v2.py - Migori (node 3) burn-in, 50 years, WITH baseline case\n'
       'management (the 2010 rows of case_management_v2.csv) and importation from day 1.\n'
       'Climate: the 10-year 2010-2019 file (3,650 days; 18,250 = 5 x 3,650 keeps the\n'
       'season aligned at the pickup start). xTLH grid 4 to 32 (13 values) x 4 seeds.\n'
       'Run from FE-2026-examples:   python run_migori_burnin_v2.py\n"""')
s = DOC + s[i1:]
rep('np.geomspace(0.5, 16.0, 13)', 'np.geomspace(4.0, 32.0, 13)')
rep('import manifest\n', 'import emodpy_malaria.interventions.treatment_seeking as cm\nimport manifest\n')
rep('SERIALIZE_YEARS = 50\n', 'SERIALIZE_YEARS = 50\nNODE = 3\n'
    'IMPORT_PER_MONTH = int(os.environ.get("IMPORT", "1"))\n'
    'INTERV_DIR = os.path.join(manifest.input_dir, "interventions")\n')
NEW = '''def fix_property_restrictions(obj):
    keys = ("Node_Property_Restrictions", "Property_Restrictions", "Property_Restrictions_Within_Node")
    if isinstance(obj, dict):
        for key, value in list(obj.items()):
            if key in keys:
                if not isinstance(value, list):
                    obj[key] = []
            else:
                fix_property_restrictions(value)
    elif isinstance(obj, list):
        for item in obj:
            fix_property_restrictions(item)


def build_camp():
    camp.set_schema(manifest.schema_file)
    cm_all = pd.read_csv(os.path.join(INTERV_DIR, "case_management_v2.csv"))
    sel = (cm_all["node_id"].astype(int) == NODE) & (cm_all["phase"] == "pickup")
    sub = cm_all[sel & (cm_all["year"] == 2010)].reset_index(drop=True)
    assert len(sub) == 4, "expected 4 case-management rows for 2010"
    targets = [{"trigger": str(sub["trigger"][r]), "coverage": float(sub["coverage"][r]),
                "agemin": float(sub["age_min"][r]), "agemax": float(sub["age_max"][r]),
                "rate": float(sub["rate"][r])} for r in range(len(sub))]
    cm.add_treatment_seeking(camp, node_ids=[NODE], start_day=1,
                             duration=SERIALIZE_YEARS * 365, drug=["Artemether", "Lumefantrine"],
                             targets=targets, broadcast_event_name="Received_Treatment")
    print("BASELINE case management (2010 rows) applied for", SERIALIZE_YEARS, "years")
    if IMPORT_PER_MONTH > 0:
        from emodpy_malaria.interventions.outbreak import add_outbreak_individual
        add_outbreak_individual(camp, start_day=1, target_num_individuals=IMPORT_PER_MONTH,
                                node_ids=[NODE], repetitions=SERIALIZE_YEARS * 365 // 30,
                                timesteps_between_repetitions=30, ignore_immunity=False)
        print("IMPORTATION:", IMPORT_PER_MONTH, "infections every 30 days")
    fix_property_restrictions(camp.campaign_dict)
    return camp
'''
rep('def build_camp():\n    camp.set_schema(manifest.schema_file)\n    return camp\n', NEW)
n = s.count('burnin_v1')
print("burnin_v1 occurrences:", n)
s = s.replace('burnin_v1', 'burnin_v2')
open("run_migori_burnin_v2.py", "w").write(s)
print("wrote run_migori_burnin_v2.py", len(s.splitlines()), "lines")
