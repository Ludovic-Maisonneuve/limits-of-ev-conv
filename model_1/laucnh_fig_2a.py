import numpy as np

eps = 0.01

tp = 1
l = 0.1
lp = 0.1
N = 10
Np = 20
d = 2 * eps
cri = 0.1 * eps
to = 0

L_cri = np.linspace(0, 1 * eps, 100)

def d_teq_tp_strong(a, b, tp, l, lp, N, Np, d, cri, s, to):

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
    G = 0.01 * np.array([[1, 0], [0, 1]])
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

L_eps_ = [1, 0.1, 0.01]

for eps_ in L_eps_:
    print('eps_', eps_, flush=True)
    a = 1 * eps_
    b = 1 * eps_
    s = 0.5 * eps * eps_
    L = [d_teq_tp_strong(a, b, tp, l, lp, N, Np, d, cri, s, to) for cri in L_cri]
    np.save('results/fig_2a/L_esp_'+str(eps_), L)


