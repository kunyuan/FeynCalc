from IO import *
import scipy.integrate as integrate
import bubble
import reduce
import matplotlib as mat
import matplotlib.pyplot as plt

plt.switch_backend('TkAgg')
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

# IsLocalField=False
IsLocalField = True

D = 3
Spin = 2

# load chemical potential shift from the file
# mu = np.loadtxt("dMu.data")
# dMu2, dMu3, dMu4 = mu[0, :]
# dMu2Err, dMu3Err, dMu4Err = mu[1, :]

# simply set all chemical potenetial shift to be zero
# dMu2, dMu3, dMu4 = -0.272, -0.109, -0.062  #rs=1, beta=40, lambda=1
# dMu2, dMu3, dMu4 = -0.261, -0.123, -0.07  #rs=1, beta=40, lambda=1.5
# dMu2, dMu3, dMu4 = -0.252, -0.129, -0.077  #rs=1, beta=40, lambda=2

dMu2, dMu3, dMu4 = -0.0936, -0.0362, -0.00  # rs=1, beta=40, lambda=1.5


Para = param(D, Spin)

DataDict, Step, Groups, ReWeight, Grids = LoadFile("data", "pid[0-9]+.dat")
KGrid = Grids["KGrid"]
print Groups

###### Calculate finite-temperature polarization ################
BubbleQ = np.zeros(len(KGrid))
for qi, q in enumerate(KGrid):
    BubbleQ[qi] = bubble.Bubble(D, Para.Beta, Spin, Para.EF, q)[0]
    # print q, BubbleQ[qi]
################################################################

Bubble = bubble.Bubble(D, Para.Beta, Spin, Para.EF, 0.0)
print "Uniform Polarization: ", Bubble[0], "+-", Bubble[1]
print "Uniform polarization for the Free electron at T=0: ", Para.Nf
Phys = Bubble[0]*len(KGrid)

EsDataDict = {}
for g in DataDict.keys():
    EsDataDict[g] = reduce.EstimateGroup(DataDict, Step, Phys, g)

for (o, key) in enumerate(sorted(EsDataDict.keys())):
    if key == (0, ):
        continue
    y = EsDataDict[key]
    print yellow("(Order: {0}, VerCT: {1}, SigmaCT: {2}) = {3:12.8f} +- {4:12.8f}".format(
        key[0], key[1], key[2], y[0][0], y[1][0]*2.0))
    EsDataDict[key] = np.array((y[0], y[1]))

Map = {}
for key in Groups:
    if key == (0, ):
        continue
    mappedkey = (key[0]+key[1], key[2])
    Map[key] = mappedkey

EsData = reduce.Reduce(EsDataDict, Map)
# print reduce.GetGroup(EsData, MappedGroups, Step, Phys, (4, 0))
# print "Mapped result: ", EsData[(4, 0)][0]

for (o, key) in enumerate(sorted(EsData.keys())):
    if key == (0, ):
        continue
    y = EsData[key]
    print green("(Order: {0}, SigmaCT: {1}) = {2:12.8f} +- {3:12.8f}".format(
        key[0], key[1], y[0][0], y[1][0]*2.0))

Each = {}
if Para.Order >= 1:
    Each[1] = EsData[(1, 0)]
if Para.Order >= 2:
    Each[2] = EsData[(2, 0)]
if Para.Order >= 3:
    Each[3] = EsData[(3, 0)]+dMu2*EsData[(1, 1)]
if Para.Order >= 4:
    Each[4] = EsData[(4, 0)]+dMu2*EsData[(2, 1)]+dMu3*EsData[(1, 1)]

Accu = {}
for o in range(1, Para.Order+1):
    Accu[o] = np.array(Each[o])
    for i in range(1, o):
        Accu[o] += Each[i]
# print Each[1].shape

fig, ax = plt.subplots()

if IsLocalField == False:
    for o in range(1, Para.Order+1):
        plt.errorbar(KGrid/Para.kF, Accu[o][0], yerr=Accu[o][1]*2.0, fmt='o-', capthick=1, capsize=4,
                     color=ColorList[o], label="Order {0}".format(o))
        # plt.errorbar(KGrid/Para.kF, Each[o][0, :], yerr=Each[o][1, :]*2.0, fmt='o-', capthick=1, capsize=4,
        #             color=ColorList[o], label="Order {0}".format(o))
        print "Order {0}: {1:12.8f} +-{2:12.8f}, Accu: {3:12.8f} +-{4:12.8f}".format(
            o, Each[o][0, 0], Each[o][1, 0]*1.0, Accu[o][0, 0], Accu[o][1, 0]*1.0)
else:
    for o in range(1, Para.Order+1):
        G = KGrid**2/8.0/np.pi*(-1.0/BubbleQ-1.0/Accu[o][0, :])

        plt.errorbar(KGrid/Para.kF, G, yerr=0.0, fmt='o-', capthick=1, capsize=4,
                     color=ColorList[o], label="Order {0}".format(o))
        # print "Order {0}: {1:12.8f} +-{2:12.8f}, Accu: {3:12.8f} +-{4:12.8f}".format(
        #     o, Each[o][0, 0], Each[o][1, 0]*1.0, Accu[o][0, 0], Accu[o][1, 0]*1.0)


# with open("chi_bubble.dat", "w") as f:
#     for ki, k in enumerate(KGrid):
#         f.write("{0} {1} {2}\n".format(k, Each[1][0, ki], Each[1][1, ki]*2.0))

ax.set_xlim([0.0, KGrid[-1]/Para.kF])
ax.set_xlabel("$q/k_F$", size=size)
ax.set_ylabel("$-P(\omega=0, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

# plt.savefig("charge_rs1_lambda1.pdf")
plt.show()
