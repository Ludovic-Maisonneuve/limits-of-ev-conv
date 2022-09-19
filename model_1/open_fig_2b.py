import matplotlib.pyplot as plt
import numpy as np

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

tp = 1
l = 0.1
lp = 0.1
N = 10
Np = 20
d = 2 * eps
cri = 0.1 * eps
to = 0
a = 1 * eps
b = 1 * eps
s = 0 * eps * eps

L_cri = np.linspace(0, 1 * eps, 100)

L_eps_1 = np.load('results/fig_2b/L_esp_1.npy')
L_eps_0_1 = np.load('results/fig_2b/L_esp_0.1.npy')
L_eps_0_01 = np.load('results/fig_2b/L_esp_0.01.npy')

def d_teq_tp(a, b, tp, l, lp, N, Np, d, cri, s, to):
    A = 2 * b * d * lp * Np / (1 + l * N + lp * Np) ** 2 + 2 * s - a * cri * Np / N

    if A <= 0:
        return float('inf')
    else:
        return 2 * s * np.abs(to - tp) / A

L_eps_approx = [d_teq_tp(a, b, tp, l, lp, N, Np, d, cri, s, to) for cri in L_cri]
for i, div in enumerate(L_eps_approx):
    if div == 0:
        j = i
y = L_cri[j]

fig, ax = plt.subplots()
plt.plot(L_cri, L_eps_1, label = r"$\tilde{\epsilon}'=1$", alpha = 1, color = 'yellowgreen')
plt.plot(L_cri, L_eps_0_1, label = r"$\tilde{\epsilon}'=0.1$", alpha = 1, color = 'mediumseagreen')
plt.plot(L_cri, L_eps_0_01, label = r"$\tilde{\epsilon}'=0.01$", alpha = 1, color = 'forestgreen')
ylim = ax.get_ylim()
plt.plot(L_cri, L_eps_approx, label = 'approx', color = 'mediumvioletred', alpha = 0.9)
plt.vlines(y, ylim[0], ylim[1], linestyles='dashed', color='mediumvioletred')
plt.ylim(ylim)
xlim = ax.get_xlim()
plt.hlines(1, xlim[0], xlim[1], linestyles='dashed', color='black', label = r'$|t_{a1}-t_{a2}|$')
plt.xlim(xlim)
plt.xlabel(r"$c_{RI}$")
plt.ylabel(r"$|\overline{t}^*_1-\overline{t}_2|$")
plt.xticks([0, 0.005, 0.01])
#plt.legend()
plt.tight_layout()
