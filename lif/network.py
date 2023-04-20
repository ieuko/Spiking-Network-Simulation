#! /usr/bin/python3

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

TAU = 20.0  # ms
V_REST = -65.0  # mV
V_RESET = -65.0  # mV
THETA = -55.0  # mV
R_M = 1.0  # MOhm
DT = 1.0  # ms
T = 1000.0  # ms; unused
NT = 1000  # (T / DT)
I_EXT = 12.0  # nA
N = 2  # of neurons
TAU_SYN = 5.0  # ms
R_SYN = 1.0  # M0hm
W = 2.0  # connection weight

v = [V_REST, V_REST - 15.]
i_syn = [0., 0.]
s = [False, False]

arr_v = np.array([])
arr_v2 = np.array([])
arr_1 = [False, False]
arr_s = np.array([[], []])

for nt in range(NT):
    t = DT * nt
    arr_v = np.append(arr_v, v[0])
    arr_v2 = np.append(arr_v2, v[1])
    for i in range(N):
        i_syn[i] = (np.exp((- DT / TAU_SYN)) * i_syn[i] + W * s[(i + 1) % 2])
    for i in range(N):
        v[i] += (DT * (- (v[i] - V_REST) + R_SYN * i_syn[i] + R_M * I_EXT) / TAU)
        s[i] = (v[i] > THETA)
        if (s[i]):
            v[i] = V_REST
            arr_1[i] = nt
    arr_s = np.append(arr_s, [[arr_1[0]], [arr_1[1]]], axis=1)


fig = plt.figure(facecolor="lightgray")

ax = fig.add_subplot(111)

x = np.arange(0, 1000, 1)

ax.set_xlim(0, 1000)
ax.set_ylim(-80, 60)

ax.plot(x, arr_v, color="blue")
ax.plot(x, arr_v2, color="silver")
ax.vlines(x=arr_s[0, :], ymin=-55.5, ymax=60, colors="blue", linewidth=0.75)
ax.vlines(x=arr_s[1, :], ymin=-55.5, ymax=60, colors="silver", linewidth=0.75)

plt.show()
