s = open("plot_county.py").read()
R = [
('CONST = 0.1635\n', 'CONST = 0.1635\nTAG = os.environ.get("TAG", "")\nLABEL = os.environ.get("LABEL", "")\nSUF = f"_{TAG}" if TAG else ""\n'),
('J = pd.read_csv(f"{COUNTY}_J_per_seed.csv");', 'J = pd.read_csv(f"{COUNTY}_J_per_seed.csv" if not TAG else f"{COUNTY}_J_{TAG}.csv");'),
('Y = pd.read_csv(f"{COUNTY}_scores_per_seed_year.csv")', 'Y = pd.read_csv(f"{COUNTY}_scores_per_seed_year.csv" if not TAG else f"{COUNTY}_scores_{TAG}.csv")'),
('fig.tight_layout(rect=(0, 0, 1, 0.96))', 'fig.tight_layout(rect=(0, 0.02, 1, 0.96))\nfig.text(0.01, 0.005, LABEL, fontsize=8.5, color=MUTE)'),
('f"{COUNTY}_calibration.{ext}"', 'f"{COUNTY}_calibration{SUF}.{ext}"'),
('f"{COUNTY}_calibration.png"))', 'f"{COUNTY}_calibration{SUF}.png"))'),
]
for a, b in R:
    assert s.count(a) == 1, a[:40]
    s = s.replace(a, b)
open("plot_county.py", "w").write(s)
print("patched")
