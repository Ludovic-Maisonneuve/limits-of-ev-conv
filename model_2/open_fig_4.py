import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import os

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

res = 100
pas = 1 / res
L_l1 = np.linspace(0 + pas, 1 - pas, res)
pas = 25 / res
L_N1 = np.linspace(0 + pas, 25 - pas, res)

blue_div = np.array([61 / 255, 139 / 255, 240 / 255])
blue = np.array([67 / 255, 114 / 255, 189 / 255])
red = np.array([159 / 255, 58 / 255, 65 / 255])
red_div = np.array([240 / 255, 58 / 255, 65 / 255])
black = np.array([0 / 255, 0 / 255, 0 / 255])
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

n = 1000
L_t = np.linspace(-1,2, n)
L = [color(t) for t in L_t]

cm = LinearSegmentedColormap.from_list(
    'cmap', L, N=n)

L_t1_impact_N1_l1_cri_0 = np.zeros((res,res))
L_t2_impact_N1_l1_cri_0 = np.zeros((res,res))
L_dis_impact_N1_l1_cri_0 = np.zeros((res,res))

L_t1_impact_N1_l1_cri_2_5 = np.zeros((res,res))
L_t2_impact_N1_l1_cri_2_5 = np.zeros((res,res))
L_dis_impact_N1_l1_cri_2_5 = np.zeros((res,res))

L_t1_impact_N1_l1_cri_5 = np.zeros((res,res))
L_t2_impact_N1_l1_cri_5 = np.zeros((res,res))
L_dis_impact_N1_l1_cri_5 = np.zeros((res,res))

thisdir = 'results/fig_4/cri0'
for r, d, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        if file[-4:] == '.npy':
            f = file.split("_")
            t1, t2 = np.load(thisdir +'/'+ file)
            L_t1_impact_N1_l1_cri_0[int(f[0]), int(f[1][:-4])] = t1
            L_t2_impact_N1_l1_cri_0[int(f[0]), int(f[1][:-4])] = t2
            L_dis_impact_N1_l1_cri_0[int(f[0]), int(f[1][:-4])] = np.abs(t2-t1)


thisdir = 'results/fig_4/cri5'
for r, d, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        if file[-4:] == '.npy':
            f = file.split("_")
            t1, t2 = np.load(thisdir +'/'+ file)
            L_t1_impact_N1_l1_cri_5[int(f[0]), int(f[1][:-4])] = t1
            L_t2_impact_N1_l1_cri_5[int(f[0]), int(f[1][:-4])] = t2
            L_dis_impact_N1_l1_cri_5[int(f[0]), int(f[1][:-4])] = np.abs(t2-t1)

gmin = -1
gmax = 2

def yc(N1):
    if l2 * N2 / N1 < 1:
        return l2 * N2 / N1
    else:
        return float('inf')

L_N1_ = np.linspace(0.001, 25, 1000)
L_yc = [yc(N1) for N1 in L_N1_]

fig, ax = plt.subplots(1)
plt.imshow(L_dis_impact_N1_l1_cri_0, extent=[0, 25, 0, 1], aspect=25, origin='lower',
           interpolation='bicubic', vmin = 0, vmax = 1.5)
plt.colorbar()
CS = ax.contour(L_N1, L_l1, L_dis_impact_N1_l1_cri_0, levels=[0.3, 0.4, 0.5, 0.6], colors='yellow')
ax.clabel(CS, inline=True, fontsize=15)
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(1)
plt.imshow(L_t1_impact_N1_l1_cri_0, extent=[0, 25, 0, 1], aspect=25, origin='lower',
           interpolation='bicubic', cmap=cm, vmin=gmin, vmax=gmax)
plt.colorbar()
ax.axvline(10, 0, 1, color='yellow', linestyle='--')
ax.axhline(0.1, 0, 25, color='yellow', linestyle='--')
ax.plot(L_N1_, L_yc, color='yellow')
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(1)
plt.imshow(L_t2_impact_N1_l1_cri_0, extent=[0, 25, 0, 1], aspect=25, origin='lower',
           interpolation='bicubic', cmap=cm, vmin=gmin, vmax=gmax)
plt.colorbar()
ax.axvline(10, 0, 1, color='yellow', linestyle='--')
ax.axhline(0.1, 0, 25, color='yellow', linestyle='--')
ax.plot(L_N1_, L_yc, color='yellow')
plt.tight_layout()
plt.show()


fig, ax = plt.subplots(1)
plt.imshow(L_dis_impact_N1_l1_cri_5, extent=[0, 25, 0, 1], aspect=25, origin='lower', vmin = 0, vmax = 1.5,
           interpolation='bicubic')
plt.colorbar()
CS = ax.contour(L_N1, L_l1, L_dis_impact_N1_l1_cri_5, levels=[0.5, 0.6, 0.7], colors='yellow')
ax.clabel(CS, inline=True, fontsize=15)
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(1)
plt.imshow(L_t1_impact_N1_l1_cri_5, extent=[0, 25, 0, 1], aspect=25, origin='lower',
           interpolation='bicubic', cmap=cm, vmin=gmin, vmax=gmax)
plt.colorbar()
ax.axvline(10, 0, 1, color='yellow', linestyle='--')
ax.axhline(0.1, 0, 25, color='yellow', linestyle='--')
ax.plot(L_N1_, L_yc, color='yellow')
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(1)
plt.imshow(L_t2_impact_N1_l1_cri_5, extent=[0, 25, 0, 1], aspect=25, origin='lower',
           interpolation='bicubic', cmap=cm, vmin=gmin, vmax=gmax)
plt.colorbar()
ax.axvline(10, 0, 1, color='yellow', linestyle='--')
ax.axhline(0.1, 0, 25, color='yellow', linestyle='--')
ax.plot(L_N1_, L_yc, color='yellow')
plt.tight_layout()
plt.show()
