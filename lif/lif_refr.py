#! /usr/bin/python3

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

TAU = 20.0  # ms
V_REST = -65.0  # mV
VRESET = -65.0  # mV
THETA = -55.0  # mV
R_M = 1.0  # MOhm
DT = 1.0  # ms
T = 1000.0  # ms; unused
NT = 1000  # (T / DT)
I_EXT = 12.0  # nA
T_REFR = 5.0  # ms; unused
NT_REFR = 5  # (T_REFR / DT)

v = V_REST
refr = 0

arr_v = np.array([])
arr_s = np.array([])

for nt in range(NT):
    t = DT * nt
    arr_v = np.append(arr_v, v)
    v += (DT * (-(v - V_REST) + R_M * I_EXT) / TAU)
    if (v > THETA):
        refr = NT_REFR
        arr_s = np.append(arr_s, nt)
    else:
        refr += -1
    if (refr > 0):
        v = V_REST


fig = plt.figure(facecolor="lightgray")

ax = fig.add_subplot(111)

x = np.arange(0, 1000, 1)

ax.set_xlim(0, 1000)
ax.set_ylim(-80, 60)

ax.plot(x, arr_v, color="blue")
ax.vlines(x=arr_s, ymin=-55.5, ymax=60, colors="blue")

plt.show()
