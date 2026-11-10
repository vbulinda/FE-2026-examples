import os, numpy as np, pandas as pd
X = {1: 0.839, 2: 0.639, 3: 0.715, 4: 0.802, 5: 0.865, 6: 0.706, 7: 0.593, 8: 0.799}
xs, taus = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0], [177, 222, 381, 639, 922, 1138]
d = os.path.join("inputs", "interventions")
t = pd.read_csv(os.path.join(d, "ITN_v5.csv"))
kill = {n: 1 / (1 + np.exp(-(0.63 + 4.0 * (x - 0.5)))) for n, x in X.items()}
tau = {n: float(np.interp(x, xs, taus)) for n, x in X.items()}
t["kill_effect"] = [f"{kill[n]:.3f}" for n in t.node_id]
t["kill_decay"] = [f"{tau[n]:.0f}" for n in t.node_id]
t.to_csv(os.path.join(d, "ITN_v6.csv"), index=False)
for n in X: print(n, X[n], round(kill[n], 3), round(tau[n]))
