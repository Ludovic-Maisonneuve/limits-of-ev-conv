import numpy as np
import os

blue_div = np.array([61 / 255, 139 / 255, 240 / 255])
blue = np.array([67 / 255, 114 / 255, 189 / 255])
red = np.array([159 / 255, 58 / 255, 65 / 255])
red_div = np.array([240 / 255, 58 / 255, 65 / 255])
gmin = -1
gmax = 2
n = 19


def color(t):
    if t == float('inf'):
        c = np.array([0 / 255, 0 / 255, 0 / 255])
        return c[0], c[1], c[2]

    i = int(n * (t + 1) / 3)

    if i < n * (-gmin) / (gmax - gmin):
        a = n * (-gmin) / (gmax - gmin)
        c = (a - i) / a * blue_div + i / a * blue
    elif i < n * (1 - gmin) / (gmax - gmin):
        a = n * (1 - gmin) / (gmax - gmin)
        b = n * (-gmin) / (gmax - gmin)
        c = (a - i) / (a - b) * blue + (i - b) / (a - b) * red
    else:
        a = n
        b = n * (1 - gmin) / (gmax - gmin)
        c = (a - i) / (a - b) * red + (i - b) / (a - b) * red_div
    return c[0], c[1], c[2]


def t1eq_t2eq(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2):
    A1 = 2 * b * d * l2 * N2 / (1 + l1 * N1 + l2 * N2) ** 2 + 2 * s1 - a1 * cri1 * N2 / N1
    A2 = 2 * b * d * l1 * N1 / (1 + l1 * N1 + l2 * N2) ** 2 + 2 * s2 - a2 * cri2 * N1 / N2

    if A1 * A2 - (2 * b * d * l2 * N2 / (1 + l1 * N1 + l2 * N2) ** 2 - a1 * cri1 * N2 / N1) * (
            2 * b * d * l1 * N1 / (1 + l1 * N1 + l2 * N2) ** 2 - a2 * cri2 * N1 / N2) > 0:
        if A1 > 0 and A2 > 0:
            C1 = (2 * b * d * l2 * N2 / (1 + l1 * N1 + l2 * N2) ** 2 - a1 * cri1 * N2 / N1) / A1
            B1 = 2 * s1 * to1 / A1
            C2 = (2 * b * d * l1 * N1 / (1 + l1 * N1 + l2 * N2) ** 2 - a2 * cri2 * N1 / N2) / A2
            B2 = 2 * s2 * to2 / A2

            t1eq = (B1 + C1 * B2) / (1 - C1 * C2)
            t2eq = (B2 + C2 * B1) / (1 - C1 * C2)
            p1eq = t1eq + cri1 * N2 / N1 * (t1eq - t2eq)
            p2eq = t2eq + cri2 * N1 / N2 * (t2eq - t1eq)

            return t1eq, t2eq

        elif A1 <= 0 and A2 <= 0:
            return float('inf'), float('inf')

        else:
            return t1_t2_simu(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2)

    else:
        return t1_t2_simu(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2)


def t1_t2_simu(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2):
    def beta(t1, p1, t2, p2):
        betat1 = - 2 * a1 * (t1 - p1) - 4 * b * d * (t1 - t2) * l2 * N2 / (
                1 + l1 * N1 + l2 * N2) ** 2 - 4 * s1 * (t1 - to1)
        betap1 = - 2 * a1 * (p1 - t1) + 2 * cri1 * a1 * (t1 - t2) * N2 / N1
        betat2 = - 2 * a2 * (t2 - p2) - 4 * b * d * (t2 - t1) * l1 * N1 / (
                1 + l1 * N1 + l2 * N2) ** 2 - 4 * s2 * (t2 - to2)
        betap2 = - 2 * a2 * (p2 - t2) + 2 * cri2 * a2 * (t2 - t1) * N1 / N2
        return np.array([betat1, betap1, betat2, betap2])

    X0 = np.array([to1, to1, to2, to2])
    X = X0
    G = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    dX = np.array([1, 1, 1, 1])
    i = 0
    while np.sqrt(dX[0] ** 2 + dX[1] ** 2 + dX[2] ** 2 + dX[3] ** 2) / 4 > 0.00000000000000001 and (
            X[0] ** 2 < 10000000 or X[1] ** 2 < 10000000 or X[2] ** 2 < 10000000 or X[
        3] ** 2 < 10000000) and i < 1000000:
        i += 1
        dX = np.dot(G, beta(X[0], X[1], X[2], X[3]))
        X = X + dX
    print(i)
    t1eq = float('inf')
    t2eq = float('inf')
    if X[0] ** 2 < 10000:
        t1eq = X[0]
    if X[2] ** 2 < 10000:
        t2eq = X[2]
    return t1eq, t2eq


eps = 0.01

a1 = 1 * eps
a2 = 1 * eps
b = 1 * eps
l1 = 0.1
l2 = 0.1
N1 = 10
N2 = 10
d = 5 * eps
cri1 = 0 * eps
cri2 = cri1
s1 = 0.5 * eps ** 2
s2 = 0.5 * eps ** 2
to1 = 0
to2 = 1

res = 250
pas = 1 / res
L_l1 = np.linspace(0 + pas, 1 - pas, res)
pas = 25 / res
L_N1 = np.linspace(0 + pas, 25 - pas, res)

print('\n cri = 0 \n', flush=True)

for i, N1 in enumerate(L_N1):
    for j, l1 in enumerate(L_l1):
        if os.path.isfile('results/fig_4/cri0/t1/' + str(j) + '_' + str(i) + '.npy') == False:
            print(i, j, flush=True)
            t1, t2 = t1eq_t2eq(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2)
            if t1 != float('inf'):
                t1 = np.min([t1, 2])
                t1 = np.max([t1, -1])
                t2 = np.min([t2, 2])
                t2 = np.max([t2, -1])
                np.save('results/fig_4/cri0/t1/' + str(j) + '_' + str(i), color(t1))
                np.save('results/fig_4/cri0/t2/' + str(j) + '_' + str(i), color(t2))

cri1 = 0.5 * eps
cri2 = cri1

print('\n cri = 0.0005 \n', flush=True)

for i, N1 in enumerate(L_N1):
    for j, l1 in enumerate(L_l1):
        if os.path.isfile('results/fig_4/cri5/t1/' + str(j) + '_' + str(i) + '.npy') == False:
            print(i, j, flush=True)

            A1 = 2 * b * d * l2 * N2 / (1 + l1 * N1 + l2 * N2) ** 2 + 2 * s1 - a1 * cri1 * N2 / N1
            A2 = 2 * b * d * l1 * N1 / (1 + l1 * N1 + l2 * N2) ** 2 + 2 * s2 - a2 * cri2 * N1 / N2

            if A1 * A2 - (2 * b * d * l2 * N2 / (1 + l1 * N1 + l2 * N2) ** 2 - a1 * cri1 * N2 / N1) * (
                    2 * b * d * l1 * N1 / (1 + l1 * N1 + l2 * N2) ** 2 - a2 * cri2 * N1 / N2) > 0 and A1 * A2 > 0:
                t1, t2 = t1eq_t2eq(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2)
                print(t1, t2)
                if t1 != float('inf'):
                    t1 = np.min([t1, 2])
                    t1 = np.max([t1, -1])
                    np.save('results/fig_4/cri5/t1/' + str(j) + '_' + str(i), color(t1))
                else:
                    np.save('results/fig_4/cri5/t1/' + str(j) + '_' + str(i), np.array([0, 0, 0]))

                if t2 != float('inf'):
                    t2 = np.min([t2, 2])
                    t2 = np.max([t2, -1])
                    np.save('results/fig_4/cri5/t2/' + str(j) + '_' + str(i), color(t2))
                else:
                    np.save('results/fig_4/cri5/t2/' + str(j) + '_' + str(i), np.array([0, 0, 0]))
