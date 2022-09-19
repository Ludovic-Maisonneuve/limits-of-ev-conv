import numpy as np
import os

blue_div = np.array([61 / 255, 139 / 255, 240 / 255])
blue = np.array([67 / 255, 114 / 255, 189 / 255])
red = np.array([159 / 255, 58 / 255, 65 / 255])
red_div = np.array([240 / 255, 58 / 255, 65 / 255])
black = np.array([0 / 255, 0 / 255, 0 / 255])
gmin = -1
gmax = 2
n = 19

def t1eq_t2eq(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2):

    return t1_t2_simu(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2)


def t1_t2_simu(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2):

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

    X0 = np.array([to1, to1, to2, to2])
    X = X0
    G = 0.01 * np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    dX = np.array([1, 1, 1, 1])
    i = 0
    while np.sqrt(dX[0] ** 2 + dX[1] ** 2 + dX[2] ** 2 + dX[3] ** 2) / 4 > 0.00000000000000001 and (
            X[0] ** 2 < 10000000 or X[1] ** 2 < 10000000 or X[2] ** 2 < 10000000 or X[
        3] ** 2 < 10000000) and i < 1000000:
        i += 1
        dX = np.dot(G, beta(X[0], X[1], X[2], X[3]))
        X = X + dX
    print(i)

    return X[0], X[2]


eps = 0.01

a1 = 1
a2 = 1
b = 1
l1 = 0.1
l2 = 0.1
N1 = 10
N2 = 10
d = 5 * eps
cri1 = 0 * eps
cri2 = cri1
s1 = 0.5 * eps
s2 = 0.5 * eps
to1 = 0
to2 = 1



res = 100
pas = 1 / res
L_l1 = np.linspace(0 + pas, 1 - pas, res)
pas = 25 / res
L_N1 = np.linspace(0 + pas, 25 - pas, res)

print('\n cri = 0 \n', flush=True)

for i, N1 in enumerate(L_N1):
    for j, l1 in enumerate(L_l1):
        if i > 1000 and os.path.isfile('results/fig_4/cri0/' + str(j) + '_' + str(i) + '.npy') == False:
            print(i, j, flush=True)
            t1, t2 = t1eq_t2eq(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2)
            np.save('results/fig_4/cri0/' + str(j) + '_' + str(i), [t1,t2])

cri1 = 0.5 * eps
cri2 = cri1

print('\n cri = 0.0005 \n', flush=True)

for i, N1 in enumerate(L_N1):
    for j, l1 in enumerate(L_l1):
        if i > 1000 and os.path.isfile('results/fig_4/cri5/' + str(j) + '_' + str(i) + '.npy') == False:
            print(i, j, flush=True)
            t1, t2 = t1eq_t2eq(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2)
            np.save('results/fig_4/cri5/' + str(j) + '_' + str(i), [t1,t2])
