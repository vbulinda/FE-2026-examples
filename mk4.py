s = open("run_bungoma_pickup_v3.py").read()
a = '"Expiration_Period_Kappa": _num(sub["discard_k"][0]),'
b = '"Expiration_Period_Lambda": _num(sub["discard_l"][0])}'
assert s.count(a) == 1 and s.count(b) == 1
s = s.replace(a, '"Expiration_Period_Kappa": float(os.environ.get("ITN_KAPPA", _num(sub["discard_k"][0]))),')
s = s.replace(b, '"Expiration_Period_Lambda": float(os.environ.get("ITN_LAMBDA", _num(sub["discard_l"][0])))}')
c = 'itn.add_scheduled_usage_dependent_bednet('
assert s.count(c) == 1
s = s.replace(c, 'print("ITN discard", year, discard)\n        ' + c)
open("run_bungoma_pickup_v4.py", "w").write(s)
print("wrote run_bungoma_pickup_v4.py", len(s.splitlines()), "lines")
