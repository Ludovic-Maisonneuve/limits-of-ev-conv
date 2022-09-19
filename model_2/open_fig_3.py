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
L_t1 = np.zeros((res, res))
L_t2 = np.zeros((res, res))

thisdir = 'results/fig_3'
for r, d, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        if file[-4:] == '.npy':
            f = file.split("_")
            t1, t2 = np.load('results/fig_3/' + f[0] + '_' + f[1][:-4] + '.npy')
            L_t1[int(f[0]), int(f[1][:-4])] = t1
            L_t2[int(f[0]), int(f[1][:-4])] = t2
            L_impact_a_b[int(f[0]), int(f[1][:-4])] = np.abs(t2-t1)

gmin = np.min(L_impact_a_b)
gmax = np.max(L_impact_a_b)

res = 100
pas = (10-0.001)/res
L_a = np.linspace(0.001+pas, 10-pas, res)
pas = (10)/res
L_b = np.linspace(0+pas, 10-pas, res)

fig, ax = plt.subplots()
plt.imshow(L_impact_a_b.T, extent=[0.001, 10, 0, 10], aspect=1, origin='lower', vmin = 0)
plt.colorbar(label=r"$|\overline{t}^*_1-\overline{t}^*_2|$")
CS = ax.contour(L_a, L_b, L_impact_a_b.T, levels=[1], colors= 'red')
plt.xlabel('female choosiness ($a$)')
plt.ylabel('predator discrimination ($b$)')
plt.tight_layout()

fig, ax = plt.subplots()
plt.imshow(L_t1.T, extent=[0.001, 10, 0, 10], aspect=1, origin='lower', vmin=gmin, vmax=gmax, cmap=cm)
plt.colorbar(label=r"$\overline{t}^*_1$")
plt.xlabel('female choosiness ($a$)')
plt.ylabel('predator discrimination ($b$)')
plt.tight_layout()

plt.figure()
plt.imshow(L_t2.T, extent=[0.001, 10, 0, 10], aspect=1, origin='lower', vmin=gmin, vmax=gmax, cmap=cm)
plt.colorbar(label=r"$\overline{t}^*_2$")
plt.xlabel('female choosiness ($a$)')
plt.ylabel('predator discrimination ($b$)')
plt.tight_layout()

