import numpy as np
import matplotlib.pyplot as plt

def formule_ts (epsilon_a, alpha, epsilon_s = 1, s0 = 1365, sigma = 5.67 * 10 ** (-8)):
    ts = (((1 - alpha) * s0 ) / (4 * sigma * (epsilon_s - epsilon_a / 2))) ** (1/4)
    return ts

def formule_ta (ts, epsilon_s = 1):
    ta = (epsilon_s * (ts**4) / 2) ** (1/4)
    return ta

#sensibilité à l'albedo (alpha)

dictionnaire_ts_albedo = {}
dictionnaire_ta_albedo = {}

for alpha in np.arange(0, 1, 0.01):
    ts = formule_ts(epsilon_a = 0.77, alpha = alpha)
    dictionnaire_ts_albedo[alpha] = ts
    ta = formule_ta(ts)
    dictionnaire_ta_albedo[alpha] = ta

plt.figure()
plt.plot(dictionnaire_ts_albedo.keys(), dictionnaire_ts_albedo.values(), label = "Température de surface")
plt.plot(dictionnaire_ta_albedo.keys(), dictionnaire_ta_albedo.values(), label = "Température de l'atmosphère")
plt.xlabel('Albédo')
plt.ylabel('Température (K)')
plt.legend()
plt.savefig("temp_albedo.png")

#sensibilité à l'emissivité de l'atmosphère (epsilon_a)

dictionnaire_ts_epsilon_a = {}
dictionnaire_ta_epsilon_a = {}

for epsilon_a in np.arange(0, 1, 0.01):
    ts = formule_ts(epsilon_a=epsilon_a, alpha=0.3)
    dictionnaire_ts_epsilon_a[epsilon_a] = ts
    ta = formule_ta(ts)
    dictionnaire_ta_epsilon_a[epsilon_a] = ta

plt.figure()
plt.plot(dictionnaire_ts_epsilon_a.keys(), dictionnaire_ts_epsilon_a.values(), label = "Température de surface")
plt.plot(dictionnaire_ta_epsilon_a.keys(), dictionnaire_ta_epsilon_a.values(), label = "Température de l'atmosphère")
plt.xlabel("Emissivité atmosphère")
plt.ylabel('Température (K)')
plt.legend()
plt.savefig("temp_ea.png")
