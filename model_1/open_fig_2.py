import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.colors import LinearSegmentedColormap

SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 15

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

res = 100

L_impact_a_b = np.zeros((res, res))

thisdir = 'results/fig_2'
for r, d, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        if file != '.DS_Store':
            f = file.split("_")
            L_impact_a_b[int(f[0]), int(f[1][:-4])] = np.load('results/fig1/' + f[0] + '_' + f[1][:-4] + '.npy')

gmin = np.min(L_impact_a_b)
gmax = np.max(L_impact_a_b)

blue_div = np.array([61 / 255, 139 / 255, 240 / 255])
blue = np.array([67 / 255, 114 / 255, 189 / 255])
red = np.array([146 / 255, 64 / 255, 79 / 255])
red = np.array([159 / 255, 58 / 255, 65 / 255])

n = 1000
L = []
for i in range(1000):
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
plt.imshow(L_impact_a_b.T, extent=[0.001, 10, 0, 10], aspect=1, origin='lower', cmap=cm, vmin=0)
plt.colorbar(label=r"$|\overline{t}^*-\overline{t}'|$")
plt.xlabel('females choosiness ($a$)')
plt.ylabel('predators discrimination ($b$)')
plt.tight_layout()

L_impact_cat = np.zeros((res, res, 3))

for i in range(res):
    for j in range(res):
        LTD = L_impact_a_b[i, j]
        if LTD == float('inf'):
            L_impact_cat[j, i, 0] = 0
            L_impact_cat[j, i, 1] = 0
            L_impact_cat[j, i, 2] = 0
        elif LTD > 1:
            L_impact_cat[j, i, 0] = 40 / 255
            L_impact_cat[j, i, 1] = 100 / 255
            L_impact_cat[j, i, 2] = 255 / 255
        elif LTD > 0.5:
            L_impact_cat[j, i, 0] = 95 / 255
            L_impact_cat[j, i, 1] = 89 / 255
            L_impact_cat[j, i, 2] = 160 / 255
        else:
            L_impact_cat[j, i, 0] = 160 / 255
            L_impact_cat[j, i, 1] = 58 / 255
            L_impact_cat[j, i, 2] = 64 / 255

fig, ax = plt.subplots(1)
ax.imshow(L_impact_cat, extent=[1, 10, 0, 10], aspect=1, origin='lower')
plt.xlabel(r'$a$')
plt.ylabel(r'$b$')
plt.tight_layout()
