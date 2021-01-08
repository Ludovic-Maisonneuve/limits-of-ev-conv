import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

SMALL_SIZE = 16
MEDIUM_SIZE = 22
BIGGER_SIZE = 25

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

eps = 0.01

## Value of cri used for fig 3
cri = 0
cri = 0.1 * eps
cri = 0.2 * eps

a = 1 * eps
b = 1 * eps
tp = 1
l = 0.1
lp = 0.1
N = 10
Np = 20
d = 5 * eps
s = 0.5 * eps ** 2
to = 0

res = 250
pas = 1 / res
L_l = np.linspace(0 + pas, 1 - pas, res)
pas = 25 / res
L_N = np.linspace(0 + pas, 25 - pas, res)


def d_teq_tp(a, b, tp, l, lp, N, Np, d, cri, s, to):
    A = 2 * b * d * lp * Np / (1 + l * N + lp * Np) ** 2 + 2 * s - a * cri * Np / N

    if A <= 0:
        return float('inf')
    else:
        return 2 * s * np.abs(to - tp) / A


L_impact_N1_l1 = np.zeros((res, res))

for i, N in enumerate(L_N):
    for j, l in enumerate(L_l):
        # print(l_N,l_lp)
        LTD = d_teq_tp(a, b, tp, l, lp, N, Np, d, cri, s, to)
        # print(LTD)
        if LTD != float('inf'):
            L_impact_N1_l1[j, i] = np.min([LTD, 2])

gmin = 0
gmax = 2

for i, N in enumerate(L_N):
    for j, l in enumerate(L_l):
        # print(l_N,l_lp)
        LTD = d_teq_tp(a, b, tp, l, lp, N, Np, d, cri, s, to)
        # print(LTD)
        if LTD == float('inf'):
            L_impact_N1_l1[j, i] = gmax

blue_div = np.array([61 / 255, 139 / 255, 240 / 255])
blue = np.array([67 / 255, 114 / 255, 189 / 255])
red = np.array([146 / 255, 64 / 255, 79 / 255])
red = np.array([159 / 255, 58 / 255, 65 / 255])

n = 11
L = []
for i in range(n):
    if i < n / gmax:
        a = int(n / gmax)
        c = (a - i) / a * red + i / a * blue
        L.append(c)
    else:
        a = int(n / gmax)
        b = int(n)
        c = (b - i) / (b - a) * blue + (i - a) / (b - a) * blue_div
        L.append(c)

cm = LinearSegmentedColormap.from_list(
    'cmap', L, N=n)

plt.figure()
plt.imshow(L_impact_N1_l1, extent=[0, 25, 0, 1], aspect=25, origin='lower', cmap=cm, vmin=0, vmax=2)
plt.colorbar(label=r"$|\overline{t}^*-\overline{t}'|$")
plt.xlabel(r'density ($N$)')
plt.ylabel(r"defence ($\lambda$)")
plt.title(r'$c_{RI} = $' + str(cri))
plt.tight_layout()
