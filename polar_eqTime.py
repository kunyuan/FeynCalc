from IO import *
import scipy.integrate as integrate
import reduce

D = 3
Spin = 2
Para = param(D, Spin)

Data, Step, Groups, ReWeight, Grids = LoadFile("data", "pid[0-9]+.dat")
assert len(Groups) == Data[0].shape[0], "group num doesn't match the data num!"
KGrid = Grids["KGrid"]
# for i in range(len(Data)):
#     print i, Data[i][1, 0]


def bubble(e):
    return Spin/8.0/np.pi**2 * \
        e**0.5/(1.0+np.cosh(Para.Beta*(e-Para.EF)))


Bubble = integrate.quad(bubble, 0.0, 100.0)
print "EqualTime Polarization: ", Bubble[0], "+-", Bubble[1]
Phys = Bubble[0]*len(KGrid)

Map = {}
for key in Groups:
    if key == (0, ):
        continue
    # if key == (2, 1, 0):
    #     continue
    mappedkey = (key[0]+key[1], key[2])
    Map[key] = mappedkey

EsData = reduce.GetData(Data, Groups, Step, Phys, Map)

fig, ax = plt.subplots()
for (o, key) in enumerate(sorted(EsData.keys())):
    y, err = EsData[key]
    print "(Order: {0}, SigamCT: {1}) = {2} +- {3}".format(
        key[0], key[1], y, err)
    EsData[key] = (np.average(y), np.average(err))

    # plt.errorbar(KGrid/Para.kF, y, yerr=err, fmt='o-', capthick=1, capsize=4,
    #              color=ColorList[o], label="Order: {0}, SigamCT: {1}".format(key[0], key[1]))

x, y = EsData[(3, 0)][0], EsData[(1, 1)][0]
xerr, yerr = EsData[(3, 0)][1], EsData[(1, 1)][1]
print yellow("Order 1 Mu CT: {0}+-{1}".format(-x/y,
                                              (abs(xerr/x)+abs(yerr/y))*abs(x/y)))


ax.set_xlim([0.0, KGrid[-1]/Para.kF])
ax.set_xlabel("$q/k_F$", size=size)
ax.set_ylabel("$-P(\\tau=0^-, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

# plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
# plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
# plt.show()
