from IO import *
import scipy.integrate as integrate

D = 3
Spin = 2
Para = param(D, Spin)

Data, Step, Groups, ReWeight, Grids = LoadFile("data", "pid[0-9]+.dat")
assert len(Groups) == Data[0].shape[0], "group num doesn't match the data num!"
KGrid = Grids["KGrid"]
print Groups
# for i in range(len(Data)):
#     print i, Data[i][1, 0]


def bubble(e):
    return Spin/8.0/np.pi**2 * \
        e**0.5/(1.0+np.cosh(Para.Beta*(e-Para.EF)))


Bubble = integrate.quad(bubble, 0.0, 100.0)
print "EqualTime Polarization: ", Bubble[0], "+-", Bubble[1]

# Phys = Bubble[0]
# Norm = [data[0, 0] for data in Data]

Phys = Bubble[0]*len(KGrid)
Norm = [sum(data[0, :]) for data in Data]
# print Phys

Accu = {}
Each = {}
for o in range(1, Para.Order+1):
    Accu[o] = []
    Each[o] = []

for (d, data) in enumerate(Data):
    polar = np.zeros([Para.Order, len(KGrid)])
    for (idx, g) in enumerate(Groups):
        if g == (0, ):
            continue
        order = g[0]+g[1]+g[2]
        if d == 1:
            print "{0}_{1} => idx {2}, with order {3}".format(
                g[0], g[1], idx, order)
        polar[order-1, :] += Data[d][idx, :]*Phys/Norm[d]

    for o in range(1, Para.Order+1):
        Each[o].append(polar[o-1, :])
        Accu[o].append(np.sum(polar[0:o, :], axis=0))

for i in range(len(Each[1])):
    print i, Each[1][i][0]

for o in Accu.keys():
    Accu[o] = Estimate(Accu[o], Norm)
    Each[o] = Estimate(Each[o], Norm)

fig, ax = plt.subplots()

for o in range(1, Para.Order+1):
    # plt.errorbar(KGrid/Para.kF, Accu[o][0], yerr=Accu[o][1]*2.0, fmt='o-', capthick=1, capsize=4,
    #              color=ColorList[o], label="Order {0}".format(o))
    plt.errorbar(KGrid/Para.kF, Each[o][0], yerr=Each[o][1]*2.0, fmt='o-', capthick=1, capsize=4,
                 color=ColorList[o], label="Order {0}".format(o))

ax.set_xlim([0.0, KGrid[-1]/Para.kF])
ax.set_xlabel("$q/k_F$", size=size)
ax.set_ylabel("$-P(\\tau=0^-, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
