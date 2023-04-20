#! /usr/bin/python3


from turtle import color
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import random

NX = 5
NY = 5
N = NX * NY

p1 = np.array([[0, 0, 1, 0, 0],
               [0, 0, 1, 0, 0],
               [1, 1, 1, 1, 1],
               [0, 0, 1, 0, 0],
               [0, 0, 1, 0, 0]])
p2 = np.array([[1, 0, 0, 0, 1],
               [0, 1, 0, 1, 0],
               [0, 0, 1, 0, 0],
               [0, 1, 0, 1, 0],
               [1, 0, 0, 0, 1]])
sort_p1 = np.array([2, 7, 10, 11, 12, 13, 14, 17, 22, 0, 1,
                   3, 4, 5, 6, 8, 9, 15, 16, 18, 19, 20, 21, 23, 24])
sort_p2 = np.array([0, 4, 6, 8, 12, 16, 18, 20, 24, 1, 2, 3,
                   5, 7, 9, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23])

TAU = 20.
W_EXC = 5.
E_REST = -65.
THETA = -55.
R_M = 16.
I_EXT = 1.5
DECAY = 0.5

NT = 1000
DT = 1.0


def wid(post, pre):
    return pre + N * post


def create_connection(w):
    for i in range(N):
        ix = int(i / NY)
        iy = i % NY
        for j in range(N):
            jx = int(j / NY)
            jy = j % NY
            if i != j:
                w[wid(i, j)] += W_EXC * p1[ix, iy] * p1[jx, jy]
                w[wid(i, j)] += W_EXC * p2[ix, iy] * p2[jx, jy]


def set_input():
    global i_ext
    i_ext = np.array([0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1,
                     1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0])


# TODO 発火頻度出力関数つくる
def print_spikes(ns):
    for x in range(NX):
        for y in range(NY):
            desc = "{0}, {1}, {2}".format(x, y, ns[y + NY * x])
            print(desc)


# TODO それぞれのニューロンの膜電位を出力する
def main():
    w = np.zeros(N * N)
    create_connection(w)
    set_input()
    v = np.zeros(N)
    for i in range(N):
        v[i] = E_REST
    g_syn = np.zeros(N)
    s = np.empty(N, dtype=np.bool_)  # spike 0 / 1
    ns = np.zeros(N)  # spike 0 / 1
    for nt in range(NT):
        for i in range(N):
            r = 0
            for j in range(N):
                x = wid(i, j)
                r += w[x] * s[j]
            g_syn[i] = DECAY * g_syn[i] + r
        for i in range(N):
            v[i] += (DT * (-(v[i] - E_REST) + g_syn[i] +
                     R_M * (i_ext[i] + random.random())) / TAU)
            s[i] = (v[i] > THETA)
            if s[i]:
                v[i] = E_REST
            ns[i] += s[i]
    print_spikes(ns)
    print(i_ext)
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    for i in range(NX):
        for j in range(NY):
            if ns[j + NY * i] < 5.0:
                color = "white"
            elif 35.0 > ns[j + NY * i] > 10.0:
                color = "lightgray"
            elif 120.0 > ns[j + NY * i] > 35.0:
                color = "gray"
            elif ns[j + NY * i] > 120.0:
                color = "black"
            else:
                color = "red"
            ax.axhspan(8 - i * 2, 10 - i * 2, j / 5, j /
                       5 + 0.2, color=color, alpha=0.5)
    plt.show()


if __name__ == "__main__":
    main()
