#! /usr/bin/python3

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

E_REST = -65.0  # mV
C = 1.0  # micro F / cm^2
G_LEAK = 0.3  # mS / cm^2
E_LEAK = 10.6 + E_REST  # mV
G_NA = 120.0  # mS / cm^2
E_NA = 115.0 + E_REST  # mV
G_K = 36.0  # mS / cm^2
E_K = -12.0 + E_REST

DT = 0.01  # 10 micros s; unused
T = 1000  # 1000 ms; unused
NT = 100000  # T / DT; unused


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


def hh(t, X, i_ext):
    v, m, h, n = X
    dmdt_r = ((1.0 / tau_m(v)) * (- m + m0(v)))
    dhdt_r = ((1.0 / tau_h(v)) * (- h + h0(v)))
    dndt_r = ((1.0 / tau_n(v)) * (- n + n0(v)))
    dvdt_ci = ((- G_LEAK * (v - E_LEAK) - G_NA * m ** 3 * h *
               (v - E_NA) - G_K * n ** 4 * (v - E_K) + i_ext) / C)
    return np.array([dvdt_ci, dmdt_r, dhdt_r, dndt_r])


t_span = [0.0, 1000.0]

x = np.arange(0, 1000, 0.01)

i_ext = 9.0  # micro A / cm^2

v = E_REST
m = m0(v)
h = h0(v)
n = n0(v)

init = [v, m, h, n]

sol = solve_ivp(hh, t_span, init, t_eval=x, args=(i_ext,))

fig = plt.figure(facecolor="lightgray")

ax = fig.add_subplot(111)

ax.set_xlim(0, 1000)
ax.set_ylim(-80, 60)

ax.plot(sol.t, sol.y[0, :], color="blue")

plt.show()
