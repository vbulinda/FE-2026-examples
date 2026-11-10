s = open("plot_county.py").read()
R = [
('top = 10 * np.ceil((top + 4) / 10)', 'top = 10 * np.ceil((top + 14) / 10)'),
('ax.legend(frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0.0, 0.88))',
 'ax.legend(frameon=False, fontsize=8.5, loc="upper left")'),
('    ax.text(2021.1, top * 0.03, "test window (not fitted)", fontsize=8.5, color=MUTE)',
 '    ax.text(2021.1, top * 0.03, "test window (not fitted)", fontsize=8.5, color=MUTE,\n            zorder=6, bbox=dict(fc="#f1f3f5", ec="none", pad=1.5))'),
('ax.text(J["xTLH"].min(), 4.5, "J = 4: close fit", color=MUTE, fontsize=8.5)',
 'ax.text(J["xTLH"].max() * 0.97, 4.6, "J = 4: close fit", color=MUTE,\n        fontsize=8.5, ha="right")'),
('ax.legend(frameon=False, fontsize=8.5, loc="upper right")',
 'ax.legend(frameon=False, fontsize=8.5, loc="upper left")'),
('    ax.text(-0.58, p + 0.8, f"{yr} survey {p}%", color=col, fontsize=8.5, ha="left")',
 '    up = yr == 2020\n    ax.text(-0.58, p + (0.8 if up else -0.8), f"{yr} survey {p}%", color=col,\n            fontsize=8.5, ha="left", va="bottom" if up else "top")'),
]
for a, b in R:
    assert s.count(a) == 1, a[:40]
    s = s.replace(a, b)
open("plot_county.py", "w").write(s)
print("patched")
