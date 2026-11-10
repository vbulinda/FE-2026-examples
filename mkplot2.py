# mkplot2.py - plot_county.py -> plot_county2.py : cleaner intervention markers and labels (Migori: LLIN lines + IRS full-effect bands in the legend; J=4 and survey labels moved into legends).
s = open("plot_county.py").read()
def rep(a, b):
    global s
    assert s.count(a) == 1, a
    s = s.replace(a, b)
rep('          "migori": [(2661, "LLIN 2017"), (2616, "IRS 2017")]}.get(COUNTY, [])',
    '          }.get(COUNTY, [])\n'
    'NETS = {"migori": [1730, 2661]}.get(COUNTY, [])   # LLIN campaigns 2014, 2017 (dotted lines)\n'
    'IRSB = {"migori": [(2616, 300), (2981, 300), (3312, 300), (3690, 150)]}.get(COUNTY, [])   # IRS rounds: start day, full-effect days (sourced-rounds runs)')
rep('    top = 10 * np.ceil((top + 14) / 10)\n', '    top = 10 * np.ceil((top + 20) / 10)\n')
rep('    ax.text(2021.1, top * 0.03, "test window (not fitted)"',
    '    for i, (d, n) in enumerate(IRSB):\n'
    '        ax.axvspan(2010 + d / 365, 2010 + (d + n) / 365, color="#9ecae1", alpha=0.3, lw=0, zorder=1,\n'
    '                   label="IRS full-effect period" if i == 0 else None)\n'
    '    for i, d in enumerate(NETS):\n'
    '        ax.axvline(2010 + d / 365, color=MUTE, ls=":", lw=1.3, label="LLIN campaign" if i == 0 else None)\n'
    '    ax.text(2021.1, top * 0.03, "test window (not fitted)"')
rep('    ax.legend(frameon=False, fontsize=8.5, loc="upper left")\n\n\nfig, axs',
    '    ax.legend(frameon=False, fontsize=8.5, loc="upper left", ncol=2)\n\n\nfig, axs')
rep('ax.axhline(4, color="#9aa5b1", lw=1, ls=":")\nax.text(J["xTLH"].max() * 0.97, 4.6, "J = 4: close fit", color=MUTE,\n        fontsize=8.5, ha="right")\n',
    'ax.axhline(4, color="#9aa5b1", lw=1, ls=":", label="J = 4: close fit")\n')
rep('    ax.axhspan(lo, hi, color=col, alpha=0.07); ax.axhline(p, color=col, ls="--", lw=1.4)\n',
    '    ax.axhspan(lo, hi, color=col, alpha=0.07); ax.axhline(p, color=col, ls="--", lw=1.4, label=f"{yr} survey {p}%")\n')
rep('    up = yr == 2020\n    ax.text(-0.58, p + (0.8 if up else -0.8), f"{yr} survey {p}%", color=col,\n            fontsize=8.5, ha="left", va="bottom" if up else "top")\n', '')
open("plot_county2.py", "w").write(s)
print("wrote plot_county2.py")
