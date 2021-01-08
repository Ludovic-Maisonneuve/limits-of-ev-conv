import matplotlib.pyplot as plt
import numpy as np

SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

eps = 0.01
L_cri = np.linspace(0, 1.5 * eps, 1000)

L_impact_cri_strong = np.load('results/fig_5_and_6/L_impact_cri_strong.npy')

fig, ax = plt.subplots()
plt.plot(L_cri, L_impact_cri_strong)
xlim = ax.get_xlim()
ylim = ax.get_ylim()
plt.xlabel(r"$c_{RI}$")
plt.ylabel(r"$|\overline{t}^*_1-\overline{t}^*_2|$")
plt.tight_layout()

L_impact_cri = np.load('results/fig_5_and_6/L_impact_cri.npy')
for i, div in enumerate(L_impact_cri):
    if div != float('inf'):
        j = i
y = L_cri[j]

fig, ax = plt.subplots()
plt.plot(L_cri, L_impact_cri)
plt.xlabel(r"$c_{RI}$")
plt.ylabel(r"$|\overline{t}^*_1-\overline{t}^*_2|$")
plt.xlim(xlim[0], xlim[1])
ylim = ax.get_ylim()
plt.vlines(y, ylim[0], ylim[1], linestyles='dashed')
plt.ylim(ylim[0], ylim[1])
ax.annotate(r'$+\infty$', xy=(0.011, 1.3), fontsize=15)
plt.tight_layout()

L_impact_cri = np.load('results/fig_5_and_6/L_impact_cri_to1_equals_to2.npy')
for i, div in enumerate(L_impact_cri):
    if div == 0:
        j = i
y = L_cri[j]

L_impact_cri_strong = np.load('results/fig_5_and_6/L_impact_cri_strong_to1_equals_to2.npy')
for i, div in enumerate(L_impact_cri_strong):
    if div > 10 ** (-5):
        j = i
        break
# y = L_cri[j]

fig, ax = plt.subplots()
plt.plot(L_cri, L_impact_cri_strong)
xlim = ax.get_xlim()
ylim = ax.get_ylim()
plt.xlabel(r"$c_{RI}$")
plt.ylabel(r"$|\overline{t}^*_1-\overline{t}^*_2|$")
plt.vlines(y, ylim[0], ylim[1], linestyles='dashed')
plt.ylim(ylim[0], ylim[1])
plt.tight_layout()

fig, ax = plt.subplots()
plt.plot(L_cri, L_impact_cri)
plt.xlabel(r"$c_{RI}$")
plt.ylabel(r"$|\overline{t}^*_1-\overline{t}^*_2|$")
plt.xlim(xlim[0], xlim[1])
plt.vlines(y, ylim[0], ylim[1], linestyles='dashed')
plt.ylim(ylim[0], ylim[1])
ax.annotate(r'$+\infty$', xy=(0.011, 0.5), fontsize=15)
plt.tight_layout()
