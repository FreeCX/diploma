from pylab import loadtxt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

plot_list = [
    ["dataB.dat", "map_B.png", "3d_B.png", "B"], 
    ["dataF1.dat", "map_F1.png", "3d_F1.png", "$|\psi_1|^2$"], 
    ["dataF2.dat", "map_F2.png", "3d_F2.png", "$|\psi_2|^2$"]
]

plt.rc('font', family="serif", size=14)
fig = plt.figure()
gcf = plt.gcf()
gcf.set_size_inches(10,8)

for load, gmap, g3d, fmt in plot_list:
# plot map
    print(">> Load file %s" % load)
    Z = loadtxt(load)
    print("-> Plot map")
    fig.clear()
    plt.xlim([400, 800])
    plt.ylim([400, 800])
    plt.xlabel('X')
    plt.ylabel('Y')
    CS = plt.contourf(Z, 256, cmap='jet')
    CS2 = plt.contour(Z, 8, colors='black', linewidth=.25)
    plt.clabel(CS2, inline=1, fontsize=10)
    fig.colorbar(CS)
    print("-> Save map to %s" % gmap)
    plt.savefig(gmap)
# plot 3d
    fig.clear()
    print("-> Plot 3d")
    plt.xlim([0,1200])
    plt.ylim([0,1200])
    X = np.arange(0, 1201, 1)
    Y = np.arange(0, 1201, 1)
    X, Y = np.meshgrid(X, Y)
    ax = fig.gca(projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel(fmt)
    surf = ax.plot_surface(X, Y, Z, rstride=6, cstride=6, cmap='jet', \
        linewidth=0, antialiased=False)
    fig.colorbar(surf)
    print("-> Save 3d to %s" % g3d)
    plt.savefig(g3d, dpi=100)

fig.clear()
print(">> Plot profile")
x1, y1 = loadtxt("profileF1.dat", unpack=True)
x2, y2 = loadtxt("profileF2.dat", unpack=True)
x3, y3 = loadtxt("profileB.dat", unpack=True)
plt.plot(x1, y1, color='black', linestyle='-.', linewidth=2, label='1')
plt.plot(x2, y2, color='black', linestyle='--', linewidth=2, label='2')
plt.plot(x3, y3, color='black', linestyle='-', linewidth=2, label='3')
plt.legend(loc = 'right')
plt.savefig("band_profile.pdf")
print("-> Save profile")