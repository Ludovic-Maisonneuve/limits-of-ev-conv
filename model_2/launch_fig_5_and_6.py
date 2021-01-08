import numpy as np

eps = 0.01

a1 = 1 * eps
a2 = 1 * eps
b = 1 * eps
l1 = 0.1
l2 = 0.1
N1 = 10
N2 = 10
d = 2 * eps
s1 = 0.5 * eps ** 2
s2 = 0.5 * eps ** 2
to1 = 1
to2 = 1

L_cri = np.linspace(0, 1.5 * eps, 1000)


def d_t1_t2(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2):
    A1 = 2 * b * d * l2 * N2 / (1 + l1 * N1 + l2 * N2) ** 2 + 2 * s1 - a1 * cri1 * N2 / N1
    A2 = 2 * b * d * l1 * N1 / (1 + l1 * N1 + l2 * N2) ** 2 + 2 * s2 - a2 * cri2 * N1 / N2

    if A1 != 0 and A2 != 0:
        C1 = (A1 - 2 * s1) / A1
        B1 = 2 * s1 * to1 / A1
        C2 = (A2 - 2 * s2) / A2
        B2 = 2 * s2 * to2 / A2

        t1eq = (B1 + C1 * B2) / (1 - C1 * C2)
        t2eq = (B2 + C2 * B1) / (1 - C1 * C2)

        Mat = np.array([[- 2 * a1 - 4 * b * d * l2 * N2 / (1 + l1 * N1 + l2 * N2) ** 2 - 2 * s1, 2 * a1,
                         4 * b * d * l2 * N2 / (1 + l1 * N1 + l2 * N2) ** 2, 0],
                        [2 * a1 + 2 * a1 * cri1 * N2 / N1, -2 * a1, - 2 * a1 * cri1 * N2 / N1, 0],
                        [4 * b * d * l1 * N1 / (1 + l1 * N1 + l2 * N2) ** 2, 0,
                         -2 * a2 - 2 * s2 - 4 * b * d * l1 * N1 / (1 + l1 * N1 + l2 * N2) ** 2, 2 * a2],
                        [-2 * a2 * cri2 * N1 / N2, 0, 2 * a2 * cri2 * N1 / N2 + 2 * a2, -2 * a2]])
        eigs = np.linalg.eig(Mat)[0]
        if eigs[0] < 0 and eigs[1] < 0 and eigs[2] < 0 and eigs[3] < 0:
            return np.abs(t1eq - t2eq)
        else:
            return float('inf')
    else:
        return float('inf')


L_impact_cri = [d_t1_t2(a1, a2, b, l1, l2, N1, N2, d, cri, cri, s1, s2, to1, to2) for cri in L_cri]
np.save('results/fig_5_and_6/L_impact_cri_to1_equals_to2', L_impact_cri)

to1 = 0
L_impact_cri = [d_t1_t2(a1, a2, b, l1, l2, N1, N2, d, cri, cri, s1, s2, to1, to2) for cri in L_cri]
np.save('results/fig_5_and_6/L_impact_cri', L_impact_cri)

a1 = 1
a2 = 1
b = 1
s1 = 0.5 * eps
s2 = 0.5 * eps
to1 = 1


def d_t1_t2(a1, a2, b, l1, l2, N1, N2, d, cri1, cri2, s1, s2, to1, to2):
    print(cri1)

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
    G = 0.01 * np.array([[1, 0.01 * a1 / (1 + 0.02 * a1), 0, 0], [0.01 * a1 / (1 + 0.02 * a1), 1, 0, 0],
                         [0, 0, 1, 0.01 * a2 / (1 + 0.02 * a2)], [0, 0, 0.01 * a2 / (1 + 0.02 * a2), 1]])
    dX = np.array([1, 1, 1, 1])
    i = 0
    while np.sqrt(dX[0] ** 2 + dX[1] ** 2 + dX[2] ** 2 + dX[3] ** 2) / 4 > 0.0000000000001 and (
            X[0] ** 2 + X[1] ** 2 + X[2] ** 2 + X[3] ** 2) / 4 < 10000000 and i < 1000000:
        i += 1
        dX = np.dot(G, beta(X[0], X[1], X[2], X[3]))
        X = X + dX
    # print(i)
    if X[0] ** 2 + X[1] ** 2 < 10000000:
        return np.abs(X[0] - X[2])
    else:
        return float('inf')


L_impact_cri_strong = [d_t1_t2(a1, a2, b, l1, l2, N1, N2, d, cri, cri, s1, s2, to1, to2) for cri in L_cri]
np.save('results/fig_5_and_6/L_impact_cri_strong_to1_equals_to2', L_impact_cri_strong)

to1 = 0
L_impact_cri_strong = [d_t1_t2(a1, a2, b, l1, l2, N1, N2, d, cri, cri, s1, s2, to1, to2) for cri in L_cri]
np.save('results/fig_5_and_6/L_impact_cri_strong', L_impact_cri_strong)
