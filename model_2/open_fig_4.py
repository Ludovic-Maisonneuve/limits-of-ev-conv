import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import os

SMALL_SIZE = 8
MEDIUM_SIZE = 16
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

eps = 0.01

a1 = 1 * eps
a2 = 1 * eps
b = 1 * eps
l1 = 0.1
l2 = 0.1
N1 = 10
N2 = 10
d = 5 * eps
cri1 = 0.5 * eps
cri2 = cri1
s1 = 0.25 * eps ** 2
s2 = 0.25 * eps ** 2
to1 = 0
to2 = 1

res = 250

L_t1_impact_N1_l1_cri_0 = np.zeros((res,res,3))
L_t2_impact_N1_l1_cri_0 = np.zeros((res,res,3))
L_t1_impact_N1_l1_cri_5 = np.zeros((res,res,3))
L_t2_impact_N1_l1_cri_5 = np.zeros((res,res,3))

thisdir = 'results/fig_4/cri0/t1'
for r, d, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        if file != '.DS_Store':
            f = file.split("_")
            L_t1_impact_N1_l1_cri_0[int(f[0]), int(f[1][:-4])] = np.load('results/fig_4/cri0/t1/' + file)

thisdir = 'results/fig_4/cri0/t2'
for r, d, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        if file != '.DS_Store':
            f = file.split("_")
            L_t2_impact_N1_l1_cri_0[int(f[0]), int(f[1][:-4])] = np.load('results/fig_4/cri0/t2/' + file)

thisdir = 'results/fig_4/cri5/t1'
for r, d, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        if file != '.DS_Store':
            f = file.split("_")
            L_t1_impact_N1_l1_cri_5[int(f[0]), int(f[1][:-4])] = np.load('results/fig_4/cri5/t1/' + file)

thisdir = 'results/fig_4/cri5/t2'
for r, d, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        if file != '.DS_Store':
            f = file.split("_")
            L_t2_impact_N1_l1_cri_5[int(f[0]), int(f[1][:-4])] = np.load('results/fig_4/cri5/t2/' + file)


gmin = -1
gmax = 2

def yc(N1):
    if l2 * N2 / N1 < 1:
        return l2 * N2 / N1
    else:
        return float('inf')

L_N1 = np.linspace(0.001, 25, 1000)
L_yc = [yc(N1) for N1 in L_N1]

fig, ax = plt.subplots(1)
plt.imshow(L_t1_impact_N1_l1_cri_0, extent=[0, 25, 0, 1], aspect=25, origin='lower', vmin=gmin, vmax=gmax,
           interpolation='bicubic')
ax.axvline(10, 0, 1, color='yellow', linestyle='--')
ax.axhline(0.1, 0, 25, color='yellow', linestyle='--')
ax.plot(L_N1, L_yc, color='yellow')
plt.tight_layout()

fig, ax = plt.subplots(1)
ax.axvline(10, 0, 1, color='yellow', linestyle='--')
ax.axhline(0.1, 0, 25, color='yellow', linestyle='--')
plt.imshow(L_t2_impact_N1_l1_cri_0, extent=[0, 25, 0, 1], aspect=25, origin='lower', vmin=gmin, vmax=gmax,
           interpolation='bicubic')
ax.plot(L_N1, L_yc, color='yellow')
plt.tight_layout()

fig, ax = plt.subplots(1)
plt.imshow(L_t1_impact_N1_l1_cri_5, extent=[0, 25, 0, 1], aspect=25, origin='lower', vmin=gmin, vmax=gmax,
           interpolation='bicubic')
ax.axvline(10, 0, 1, color='yellow', linestyle='--')
ax.axhline(0.1, 0, 25, color='yellow', linestyle='--')
ax.plot(L_N1, L_yc, color='yellow')
plt.tight_layout()

fig, ax = plt.subplots(1)
ax.axvline(10, 0, 1, color='yellow', linestyle='--')
ax.axhline(0.1, 0, 25, color='yellow', linestyle='--')
plt.imshow(L_t2_impact_N1_l1_cri_5, extent=[0, 25, 0, 1], aspect=25, origin='lower', vmin=gmin, vmax=gmax,
           interpolation='bicubic')
ax.plot(L_N1, L_yc, color='yellow')
plt.tight_layout()

blue_div = np.array([61 / 255, 139 / 255, 240 / 255])
blue = np.array([67 / 255, 114 / 255, 189 / 255])
red = np.array([159 / 255, 58 / 255, 65 / 255])
red_div = np.array([240 / 255, 58 / 255, 65 / 255])
gmin = -1
gmax = 2
n = 19

def color(t):
    i = int(n * (t + 1) / 3)
    if t == float('inf'):
        c = np.array([0 / 255, 0 / 255, 0 / 255])
    elif i < n * (-gmin) / (gmax - gmin):
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


L = []
for t in np.linspace(-1, 2, n):
    L.append(color(t))

cm = LinearSegmentedColormap.from_list(
    'cmap', L, N=n)

plt.figure()
plt.imshow(np.zeros((res, res)).T, extent=[0.001, 10, 0, 10], aspect=1, origin='lower', cmap=cm, vmin=-1, vmax=2)
plt.colorbar()
plt.tight_layout()
