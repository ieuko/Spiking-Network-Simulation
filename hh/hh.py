#! /usr/bin/python3
#! This program takes a little time to process

import numpy as np
import matplotlib.pyplot as plt

E_REST = -65.0  # mV
C = 1.0  # micro F / cm^2
G_LEAK = 0.3  # mS / cm^2
E_LEAK = 10.6 + E_REST  # mV
G_NA = 120.0  # mS / cm^2
E_NA = 115.0 + E_REST  # mV
G_K = 36.0  # mS / cm^2
E_K = -12.0 + E_REST

DT = 0.01  # 10 micros s
T = 1000  # 1000 ms; unused
NT = 100000  # T / DT


def alpha_m(v):
    alpha_m_k = ((2.5 - 0.1 * (v - E_REST)) /
                 ((np.e ** (2.5 - 0.1 * (v - E_REST))) - 1.0))
    return alpha_m_k


def beta_m(v):
    beta_m_k = (4.0 * np.e ** (-(v - E_REST) / 18.0))
    return beta_m_k


def alpha_h(v):
    alpha_h_k = (0.07 * np.e ** (-(v - E_REST) / 20.0))
    return alpha_h_k


def beta_h(v):
    beta_h_k = (1.0 / (np.e ** (3.0 - 0.1 * (v - E_REST)) + 1.0))
    return beta_h_k


def alpha_n(v):
    alpha_n_k = ((0.1 - 0.01 * (v - E_REST)) /
                 (np.e ** (1 - 0.1 * (v - E_REST)) - 1.0))
    return alpha_n_k


def beta_n(v):
    beta_n_k = (0.125 * np.e ** (-(v - E_REST) / 80.0))
    return beta_n_k


def m0(v):
    m0_r = (alpha_m(v) / (alpha_m(v) + beta_m(v)))
    return m0_r


def h0(v):
    h0_r = (alpha_h(v) / (alpha_h(v) + beta_h(v)))
    return h0_r


def n0(v):
    n0_r = (alpha_n(v) / (alpha_n(v) + beta_n(v)))
    return n0_r


def tau_m(v):
    tau_m_t = (1. / (alpha_m(v) + beta_m(v)))
    return tau_m_t


def tau_h(v):
    tau_h_t = (1. / (alpha_h(v) + beta_h(v)))
    return tau_h_t


def tau_n(v):
    tau_n_t = (1. / (alpha_n(v) + beta_n(v)))
    return tau_n_t


def dmdt(v, m):
    dmdt_r = ((1.0 / tau_m(v)) * (- m + m0(v)))
    return dmdt_r


def dhdt(v, h):
    dhdt_r = ((1.0 / tau_h(v)) * (- h + h0(v)))
    return dhdt_r


def dndt(v, n):
    dndt_r = ((1.0 / tau_n(v)) * (- n + n0(v)))
    return dndt_r


def dvdt(v, m, h, n, i_ext):
    dvdt_ci = ((- G_LEAK * (v - E_LEAK) - G_NA * m ** 3 * h *
               (v - E_NA) - G_K * n ** 4 * (v - E_K) + i_ext) / C)
    return dvdt_ci


v = E_REST
m = m0(v)
h = h0(v)
n = n0(v)

i_ext = 9.0  # micro A / cm^2

v_arr = np.array([])

for i in range(NT):
    dmdt1 = dmdt(v, m)
    dhdt1 = dhdt(v, h)
    dndt1 = dndt(v, n)
    dvdt1 = dvdt(v, m, h, n, i_ext)
    dmdt2 = dmdt(v + .5 * DT * dvdt1, m + .5 * DT * dmdt1)
    dhdt2 = dhdt(v + .5 * DT * dvdt1, h + .5 * DT * dhdt1)
    dndt2 = dndt(v + .5 * DT * dvdt1, n + .5 * DT * dndt1)
    dvdt2 = dvdt(v + .5 * DT * dvdt1, m + .5 * DT * dmdt1,
                 h + .5 * DT * dhdt1, n + .5 * DT * dndt1, i_ext)
    dmdt3 = dmdt(v + .5 * DT * dvdt2, m + .5 * DT * dmdt2)
    dhdt3 = dhdt(v + .5 * DT * dvdt2, h + .5 * DT * dhdt2)
    dndt3 = dndt(v + .5 * DT * dvdt2, n + .5 * DT * dndt2)
    dvdt3 = dvdt(v + .5 * DT * dvdt2, m + .5 * DT * dmdt2,
                 h + .5 * DT * dhdt2, n + .5 * DT * dndt2, i_ext)
    dmdt4 = dmdt(v + .5 * DT * dvdt3, m + .5 * DT * dmdt3)
    dhdt4 = dhdt(v + .5 * DT * dvdt3, h + .5 * DT * dhdt3)
    dndt4 = dndt(v + .5 * DT * dvdt3, n + .5 * DT * dndt3)
    dvdt4 = dvdt(v + .5 * DT * dvdt3, m + .5 * DT * dmdt3,
                 h + .5 * DT * dhdt3, n + .5 * DT * dndt3, i_ext)
    m += DT * (dmdt1 + 2 * dmdt2 + 2 * dmdt3 + dmdt4) / 6.
    h += DT * (dhdt1 + 2 * dhdt2 + 2 * dhdt3 + dhdt4) / 6.
    n += DT * (dndt1 + 2 * dndt2 + 2 * dndt3 + dndt4) / 6.
    v += DT * (dvdt1 + 2 * dvdt2 + 2 * dvdt3 + dvdt4) / 6.
    v_arr = np.append(v_arr, v)


fig = plt.figure(facecolor="lightgray")

ax = fig.add_subplot(111)

x = np.arange(0, 1000, 0.01)

ax.set_xlim(0, 100)
ax.set_ylim(-80, 60)

ax.plot(x, v_arr, color="blue")

plt.show()
