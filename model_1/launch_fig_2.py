import numpy as np
import os

eps = 0.01

a = 1
b = 1
tp = 1
l = 0.1
lp = 0.1
N = 10
Np = 20
d = 2 * eps
cri = 0.2 * eps
s = 0.5 * eps
to = 0

res = 100
pas = (10-0.001)/res
L_a = np.linspace(0.001+pas, 10-pas, res)
pas = (10)/res
L_b = np.linspace(0+pas, 10-pas, res)

def d_teq_tp(a, b, tp, l, lp, N, Np, d, cri, s, to):

    def beta(t, p):
        betat = - 2 * a * (t - p) - 4 * b * (t - tp) * lp * Np * np.exp(- b * (t - tp) ** 2) * (
                1 / (1 + l * N + lp * Np * np.exp(- b * (t - tp) ** 2) - d) - 1 / (
                1 + l * N + lp * Np * np.exp(- b * (t - tp) ** 2))) - 4 * s * (t - to)
        betap = - 2 * a * (p - t) + 2 * cri * a * (t - tp) / (cri + np.exp(a * (t - tp) * (2 * p - t - tp)) * N / Np)
        return np.array([betat, betap])

    teq = to
    peq = to

    X0 = np.array([teq, peq])
    X = X0
    G = 0.01 * np.array([[1, 0.01 * a / (1 + 0.02 * a)], [0.01 * a / (1 + 0.02 * a), 1]])
    dX = np.array([1, 1])
    i = 0
    while np.sqrt(dX[0]**2 + dX[1]**2) > 0.000000000000001 and X[0]**2 + X[1]**2 < 10000000:
        i+= 1
        dX = np.dot(G, beta(X[0], X[1]))
        X = X + dX
    print(i)
    if X[0]**2 + X[1]**2 < 10000000:
        return np.abs(X[0] - tp)
    else:
        return float('inf')

L_impact_a_b = np.zeros((len(L_a),len(L_b)))

for i, a_ in enumerate(L_a):
    for j, b_ in enumerate(L_b):
        if os.path.isfile('results/fig_2/' + str(i) + '_' + str(j) + '.npy') == False:
            print(i,j,flush=True)
            div = d_teq_tp(a_, b_, tp, l, lp, N, Np, d, cri, s, to)
            np.save('results/fig_2/' + str(i) + '_' + str(j),div)
