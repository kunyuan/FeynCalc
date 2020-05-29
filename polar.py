from IO import *
import scipy.integrate as integrate

D = 3
Spin = 2
Para = param(D, Spin)

Data, Groups, Grids = LoadFile("data", "pid[0-9]+.dat")
KGrid = Grids["KGrid"]
Bubble = 0.0971613  # 3D, T=0.04Ef, rs=1


def bubble(e):
    return Para.Beta*Spin/8.0/np.pi**2 * \
        e**0.5/(1.0+np.cosh(Para.Beta*(e-Para.EF)))


Bubble = integrate.quad(bubble, 0.0, 100.0)
print "Static Polarization: ", Bubble[0], "+-", Bubble[1]

Phys = Bubble[0]*len(KGrid)
Norm = [sum(data[0, :]) for data in Data]

Accu = {}
for o in range(1, Para.Order+1):
    Accu[o] = [np.zeros(len(KGrid)), ]*len(Data)

for (idx, g) in enumerate(Groups):
    if g == (0, ):
        continue
    order = g[0]+g[1]
    print "{0}_{1} => {2}".format(g[0], g[1], order)
    for (di, d) in enumerate(Data):
        # print Accu[o][di].shape
        # print d[order, :].shape
        # print Norm[di]
        Accu[o][di] += d[order, :]/Norm[di]*Phys

Polar = {}
for o in Accu.keys():
    Polar[o] = Estimate(Accu[o], Norm)

print Polar[1]
print Polar[2]

fig, ax = plt.subplots()


# def ErrorPlot(p, d, color, marker, label=None, size=4, shift=False):
#     data = np.array(d)
#     data[:, 0] /= kF
#     if shift:
#         data[:, 1] -= data[0, 1]
#     p.plot(data[:, 0], data[:, 1], marker=marker, c=color, label=label,
#            lw=1, markeredgecolor="None", linestyle="--", markersize=size)
# p.errorbar(data[:,0],data[:,1], yerr=data[:,2], c=color, ecolor=color, capsize=0, linestyle="None")
# p.fill_between(data[:,0], data[:,1]-data[:,2], data[:,1]+data[:,2], alpha=0.5, facecolor=color, edgecolor=color)


# for i in range(0, len(ScanOrder)):
#     o = ScanOrder[i]
#     ErrorPlot(ax, DataOrderByOrder[o],
#               ColorList[i], 's', "Order {0}".format(o))
# ErrorPlot(ax, DataAtOrder[o], ColorList[i], 's', "Order {0}".format(o))

# ErrorPlot(ax, Data[1][0], 'k', 's', "Diag 1")
# ErrorPlot(ax, tmp, 'm', 's', "Diag 3+c 1")
# ErrorPlot(ax, DataAll[3], 'k', 'o', "Order 3")

# ErrorPlot(ax, Data[2][1], 'g', 'o', "Order 3 counterbubble 1")
# ErrorPlot(ax, Data[2][2], 'g', '*', "Order 3 counterbubble 2")
# ErrorPlot(ax, Data[2][3], 'g', '>', "Order 3 counterbubble 3")

# ErrorPlot(ax, Data[1][1], 'olive', 'o', "Order 3 shift 1")
# ErrorPlot(ax, Data[1][2], 'olive', '*', "Order 3 shift 2")

# ErrorPlot(ax, Data[3][0], 'k', 's', "Diag 1", shift=True)
# ErrorPlot(ax, Data[3][1], 'g', 's', "Diag 2", shift=True)
# ErrorPlot(ax, Data[3][2], 'r', '*', "Diag 3", shift=True)
# ErrorPlot(ax, Data[3][3], 'b', 's', "Diag 4", shift=True)
# ErrorPlot(ax, Data[3][4], 'olive', '*', "Diag 5", shift=True)
# ErrorPlot(ax, Data[3][5], 'm', 's', "Diag 6", shift=True)
# ErrorPlot(ax, Data[3][6], 'c', '*', "Diag 7", shift=True)

# ErrorPlot(ax, Data[5], 'g', 's', "Diag 6")


# x = np.arange(0, 0.2, 0.001)
# y = 0.5*x**w
# ax.plot(x,y,'k-', lw=2)

ax.set_xlim([0.0, KGrid[-1]/Para.kF])
# ax.set_xticks([0.0,0.04,0.08,0.12])
# ax.set_yticks([0.35,0.4,0.45,0.5])
# ax.set_ylim([-0.02, 0.125])
# ax.set_ylim([0.07, 0.125])
ax.set_xlabel("$q/k_F$", size=size)
# ax.xaxis.set_label_coords(0.97, -0.01)
# # ax.yaxis.set_label_coords(0.97, -0.01)
# ax.text(-0.012,0.52, "$-I$", fontsize=size)
ax.set_ylabel("$-P(\omega=0, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
