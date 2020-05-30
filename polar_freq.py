from IO import *
import scipy.integrate as integrate
import reduce

D = 3
Spin = 2
Para = param(D, Spin)

Data, Step, Groups, ReWeight, Grids = LoadFile("data", "pid[0-9]+.dat")
assert len(Groups) == Data[0].shape[0], "group num doesn't match the data num!"
KGrid = Grids["KGrid"]
print Groups


def bubble(e):
    return Para.Beta*Spin/8.0/np.pi**2 * \
        e**0.5/(1.0+np.cosh(Para.Beta*(e-Para.EF)))


Bubble = integrate.quad(bubble, 0.0, 100.0)
print "Static Polarization: ", Bubble[0], "+-", Bubble[1]

# Phys = Bubble[0]
# Norm = [data[0, 0] for data in Data]

Phys = Bubble[0]*len(KGrid)
Norm = [sum(data[0, :]) for data in Data]
# print Phys

Map = {}
for key in Groups:
    if key == (0, ):
        continue
    mappedkey = (key[0]+key[1], key[2])
    Map[key] = mappedkey

# del Map[(3, 0, 0)]

EsData = reduce.GetData(Data, Groups, Step, Phys, Map)

# dMu2, dMu3, dMu4 = -0.0, 0.0, 0.0
dMu2, dMu3, dMu4 = -0.2729, 0.0, 0.0
# dMu2, dMu3, dMu4 = -0.267, 0.0, 0.0
dErr2, dErr3, dErr4 = 0.001, 0.0, 0.0
Each = {}
Each[1] = EsData[(1, 0)]
Each[2] = EsData[(2, 0)]
Each[3] = EsData[(3, 0)]+dMu2*EsData[(1, 1)]

Accu = {}
for o in range(1, Para.Order+1):
    Accu[o] = np.array(Each[o])
    for i in range(1, o):
        Accu[o] += Each[i]

fig, ax = plt.subplots()

for o in range(1, Para.Order+1):
    # plt.errorbar(KGrid/Para.kF, Accu[o][0], yerr=Accu[o][1]*2.0, fmt='o-', capthick=1, capsize=4,
    #              color=ColorList[o], label="Order {0}".format(o))
    plt.errorbar(KGrid/Para.kF, Each[o][0, :], yerr=Each[o][1, :]*2.0, fmt='o-', capthick=1, capsize=4,
                 color=ColorList[o], label="Order {0}".format(o))
    print "Order {0}: {1}+-{2}, Accu: {3}+-{4}".format(
        o, Each[o][0, 0], Each[o][1, 0]*1.0, Accu[o][0, 0], Accu[o][1, 0]*1.0)

ax.set_xlim([0.0, KGrid[-1]/Para.kF])
ax.set_xlabel("$q/k_F$", size=size)
ax.set_ylabel("$-P(\omega=0, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
