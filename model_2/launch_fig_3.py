import numpy as np
import os

eps = 0.01

#a1 = 1 * eps
#a2 = 1 * eps
#b = 1 * eps
l1 = 0.1
l2 = 0.1
N1 = 20
N2 = 20
d = 2 * eps
cri1 = 0.2 * eps
cri2 = cri1
s1 = 0.5 * eps
s2 = 0.5 * eps
to1 = 0
to2 = 1

res = 100
pas = (10-0.001)/res
L_a = np.linspace(0.001+pas, 10-pas, res)
pas = (10)/res
L_b = np.linspace(0+pas, 10-pas, res)

def t1_t2(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2):

    def beta(t1, p1, t2, p2):
        betat1 = - 2 * a1 * (t1 - p1) - 4 * b * (t1 - t2) * l2 * N2 * np.exp(- b * (t1 - t2) ** 2) * (
                1 / (1 + l1 * N1 + l2 * N2 * np.exp(- b * (t1 - t2) ** 2) - d) - 1 / (
                1 + l1 * N1 + l2 * N2 * np.exp(- b * (t1 - t2) ** 2))) - 4 * s1 * (t1 - to1)
        betap1 = - 2 * a1 * (p1 - t1) + 2 * cri1 * a1 * (t1 - t2) / (
                cri1 + np.exp(a1 * (t1 - t2) * (2 * p1 - t1 - t2)) * N1 / N2)
        betat2 = - 2 * a2 * (t2 - p2) - 4 * b * (t2 - t1) * l1 * N1 * np.exp(- b * (t2 - t1) ** 2) * (
                1 / (1 + l1 * N1 * np.exp(- b * (t1 - t2) ** 2) + l2 * N2 - d) - 1 / (
                1 + l1 * N1 * np.exp(- b * (t1 - t2) ** 2) + l2 * N2)) - 4 * s2 * (t2 - to2)
        betap2 = - 2 * a2 * (p2 - t2) + 2 * cri2 * a2 * (t2 - t1) / (
                cri2 + np.exp(a2 * (t2 - t1) * (2 * p2 - t1 - t2)) * N2 / N1)
        return np.array([betat1, betap1, betat2, betap2])

    X0 = np.array([to1, to1, to2 + eps, to2 + eps])
    X = X0
    G = 0.01 * np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                         [0, 0, 1, 0], [0, 0, 0, 1]])
    dX = np.array([1, 1, 1, 1])
    i = 0
    while np.sqrt(dX[0] ** 2 + dX[1] ** 2 + dX[2] ** 2 + dX[3] ** 2) / 4 > 0.0000000000001 and (
            X[0] ** 2 + X[1] ** 2 + X[2] ** 2 + X[3] ** 2) / 4 < 10000000 and i < 1000000:
        i += 1
        dX = np.dot(G, beta(X[0], X[1], X[2], X[3]))
        X = X + dX
    # print(i)
    # if X[0] ** 2 + X[1] ** 2 < 10000000:
    #     return np.abs(X[0] - X[2])
    # else:
    #     return float('inf')
    return X[0], X[2]

L_impact_a_b = np.zeros((len(L_a),len(L_b)))

for i, a in enumerate(L_a):
    for j, b in enumerate(L_b):
        if i > -1 and os.path.isfile('results/fig_3/' + str(i) + '_' + str(j) + '.npy') == False:
            print(i,j,flush=True)
            X = t1_t2(a, a, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2)
            np.save('results/fig_3/' + str(i) + '_' + str(j), X)