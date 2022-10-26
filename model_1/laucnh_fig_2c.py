import numpy as np

eps = 0.01

tp = 1
l = 0.1
lp = 0.1
N = 10
Np = 20
d = 2 * eps
cri = 0.1 * eps
to = 1

L_cri = np.linspace(0, 1 * eps, 100)

def d_teq_tp_strong(a, b, tp, l, lp, N, Np, d, cri, s, to): #function returning the phenotypic distance between the two species at equilibrium

    def beta(t, p): #function computing the selection vector
        betat = - 2 * a * (t - p) - 4 * b * (t - tp) * lp * Np * np.exp(- b * (t - tp) ** 2) * (
                1 / (1 + l * N + lp * Np * np.exp(- b * (t - tp) ** 2) - d) - 1 / (
                1 + l * N + lp * Np * np.exp(- b * (t - tp) ** 2))) - 4 * s * (t - to)
        betap = - 2 * a * (p - t) + 2 * cri * a * (t - tp) / (cri + np.exp(a * (t - tp) * (2 * p - t - tp)) * N / Np)
        return np.array([betat, betap])

    teq = to+0.01 #set initial trait in species 1
    peq = to+0.01 #set initial preference in species 1

    X0 = np.array([teq, peq]) #vector with initial trait and preference
    X = X0
    G = 0.01 * np.array([[1, 0], [0, 1]]) #define matrix of genetic covariance
    dX = np.array([1, 1])
    i = 0
    while np.sqrt(dX[0]**2 + dX[1]**2) > 0.000000000000001 and X[0]**2 + X[1]**2 < 10000000 and i < 1000000: #while equilibrium is not reached
        i+= 1
        dX = np.dot(G, beta(X[0], X[1]))
        X = X + dX #evolutionary dynamics of the mean trait and preference values
    print(i, flush=True)
    if X[0]**2 + X[1]**2 < 10000000:
        return np.abs(X[0] - tp)
    else:
        return float('inf')

def d_teq_tp_strong_10(a, b, tp, l, lp, N, Np, d, cri, s, to):

    def beta(t, p):
        betat = - 2 * a * (t - p) - 4 * b * (t - tp) * lp * Np * np.exp(- b * (t - tp) ** 2) * (
                1 / (1 + l * N + lp * Np * np.exp(- b * (t - tp) ** 2) - d) - 1 / (
                1 + l * N + lp * Np * np.exp(- b * (t - tp) ** 2))) - 4 * s * (t - to)
        betap = - 2 * a * (p - t) + 2 * cri * a * (t - tp) / (cri + np.exp(a * (t - tp) * (2 * p - t - tp)) * N / Np)
        return np.array([betat, betap])

    teq = to+0.01
    peq = to+0.01

    X0 = np.array([teq, peq])
    X = X0
    G = 0.01 * np.array([[1, 0], [0, 1]])
    dX = np.array([1, 1])
    i = 0
    while np.sqrt(dX[0]**2 + dX[1]**2) > 0.000000000000001 and X[0]**2 + X[1]**2 < 10000000 and i < 1000000:
        i+= 1
        dX = np.dot(G, beta(X[0], X[1]))
        X = X + 10 * dX
    print(i, flush=True)
    if X[0]**2 + X[1]**2 < 10000000:
        return np.abs(X[0] - tp)
    else:
        return float('inf')

def d_teq_tp_strong_100(a, b, tp, l, lp, N, Np, d, cri, s, to):

    def beta(t, p):
        betat = - 2 * a * (t - p) - 4 * b * (t - tp) * lp * Np * np.exp(- b * (t - tp) ** 2) * (
                1 / (1 + l * N + lp * Np * np.exp(- b * (t - tp) ** 2) - d) - 1 / (
                1 + l * N + lp * Np * np.exp(- b * (t - tp) ** 2))) - 4 * s * (t - to)
        betap = - 2 * a * (p - t) + 2 * cri * a * (t - tp) / (cri + np.exp(a * (t - tp) * (2 * p - t - tp)) * N / Np)
        return np.array([betat, betap])

    teq = to+0.01
    peq = to+0.01

    X0 = np.array([teq, peq])
    X = X0
    G = 0.01 * np.array([[1, 0], [0, 1]])
    dX = np.array([1, 1])
    i = 0
    while np.sqrt(dX[0]**2 + dX[1]**2) > 0.000000000000001 and X[0]**2 + X[1]**2 < 10000000 and i < 1000000:
        i+= 1
        dX = np.dot(G, beta(X[0], X[1]))
        X = X + 100 * dX
    print(i, flush=True)
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
    if eps_ == 1:
        L = [d_teq_tp_strong(a, b, tp, l, lp, N, Np, d, cri, s, to) for cri in L_cri]
    elif eps_ == 0.1:
        L = [d_teq_tp_strong_10(a, b, tp, l, lp, N, Np, d, cri, s, to) for cri in L_cri]
    else:
        L = [d_teq_tp_strong_100(a, b, tp, l, lp, N, Np, d, cri, s, to) for cri in L_cri]
    np.save('results/fig_2c/L_esp_'+str(eps_), L)
