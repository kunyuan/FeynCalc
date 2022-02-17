# Calculate dMu in counterterms
#!/usr/bin/python
from typing import OrderedDict
from IO import *
import bubble
import scipy.integrate as integrate
import reduce
import matplotlib as mat
#import matplotlib.pyplot as plt

# plt.switch_backend('TkAgg')
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

D = 3
Spin = 2

inlist = open("./inlist", "r")
for index, eachline in enumerate(inlist):
    para = eachline.split()
    if len(para) == 0:
        break
    with open("./parameter", "w") as file:
        parameters = ' '.join(para[:-2])
        file.write(parameters+"\n\n")
        file.write(
            "#Order, Beta, rs, Mass2, Lambda, MinExtMom(*kF), MaxExtMom, TotalStep(*1e6)")

    Para = param(D, Spin)

    DataDict, Step, Groups, ReWeight, Grids = LoadFile(
        Para.DataFolder, "pid[0-9]+.dat")
    KGrid = Grids["KGrid"]
    # for i in range(len(Data)):
    #     print i, Data[i][1, 0]

    Bubble = bubble.Bubble(D, Para.Beta, Spin, Para.kF, 0.0)
    print("Uniform Polarization: ", Bubble[0], "+-", Bubble[1])
    print("Uniform polarization for the Free electron at T=0: ", Para.Nf)
    Phys = Bubble[0]*len(KGrid)
    # Phys = Para.Nf*len(KGrid)

    EsDataDict = {}
    for g in DataDict.keys():
        EsDataDict[g] = reduce.EstimateGroup(DataDict, Step, Phys, g)

    for (o, key) in enumerate(sorted(EsDataDict.keys())):
        if key == (0, ):
            continue
        y = EsDataDict[key]
        print(yellow("(Order: {0}, VerCT: {1}, SigamCT: {2}) = {3:12.6f} +- {4:12.6f}".format(
            key[0], key[1], key[2], np.average(y[0]), np.average(y[1])*2.0)))
        EsDataDict[key] = np.array((np.average(y[0]), np.average(y[1])))

    # print EsDataDict[(4, 0, 0)][0]
    # print EsDataDict[(3, 1, 0)][0]
    # print EsDataDict[(2, 2, 0)][0]

    Map = {}
    for key in Groups:
        if key == (0, ):
            continue
        # if key == (2, 1, 0):
        #     continue
        mappedkey = (key[0]+key[1], key[2])
        Map[key] = mappedkey

    # EsData, MappedGroups = reduce.GetData(Data, Groups, Step, Phys, Map)
    EsData = reduce.Reduce(EsDataDict, Map)
    # print reduce.GetGroup(EsData, MappedGroups, Step, Phys, (4, 0))
    print("Mapped result: ", EsData[(4, 0)][0])

    #fig, ax = plt.subplots()
    for (o, key) in enumerate(sorted(EsData.keys())):
        if key == (0, ):
            continue
        y = EsData[key]
        print(green("(Order: {0}, SigamCT: {1}) = {2:12.6f} +- {3:12.6f}".format(
            key[0], key[1], np.average(y[0]), np.average(y[1])*2.0)))

    # plt.errorbar(KGrid/Para.kF, y, yerr=err, fmt='o-', capthick=1, capsize=4,
    #              color=ColorList[o], label="Order: {0}, SigamCT: {1}".format(key[0], key[1]))

    s30, s11 = EsData[(3, 0)][0], EsData[(1, 1)][0]
    e30, e11 = EsData[(3, 0)][1], EsData[(1, 1)][1]
    # print('EsData(1,1)', s11, 'EsData(3,0)', s30)
    dMu2, dMu2Err = -s30/s11, (abs(e30/s30)+abs(e11/s11))*abs(s30/s11)*2.0
    print(yellow("Order 2 Mu CT: {0}+-{1}".format(dMu2, dMu2Err)))

    s40, s21 = EsData[(4, 0)][0], EsData[(2, 1)][0]
    e40, e21 = EsData[(4, 0)][1], EsData[(2, 1)][1]
    # print('EsData(2,1)', s21, 'EsData(4,0)', s40)
    dMu3 = -(s40+s21*dMu2)/s11
    dMu3Err = (abs((e40+e21*dMu2)/(s40+s21*dMu2))+abs(e11/s11))*abs(dMu3)*2.0
    print(yellow("Order 3 Mu CT: {0}+-{1}".format(dMu3, dMu3Err)))

    s50, s12, s31 = EsData[(5, 0)][0], EsData[(1, 2)][0], EsData[(3, 1)][0]
    e50, e12, e31 = EsData[(5, 0)][1], EsData[(1, 2)][1], EsData[(3, 1)][1]
    # print('EsData(1,2)', s12, 'EsData(3,1)', s31, 'EsData(5,0)', s50)
    dMu4 = -(s50+s12*dMu2**2+s21*dMu3+s31*dMu2)/s11
    dMu4Err = (abs((e50+e12*dMu2**2+e21*dMu3+e31*dMu2) /
                   (s50+s12*dMu2**2+s21*dMu3+s31*dMu2))+abs(e11/s11))*abs(dMu4)*2.0
    print(yellow("Order 4 Mu CT: {0}+-{1}".format(dMu4, dMu4Err)))

    with open(Para.DataFolder+"/parameter", "r") as file:
        para = file.readline().split(" ")
        beta = float(para[1])
        rs = float(para[2])
        lam = float(para[4])
        order = int(para[0])

    filename = "dMu/dMu_beta{0}_rs{1}_lam{2}".format(beta, rs, lam, order)
    with open(filename, "w") as f:
        print('Save dMu.data in '+filename)
        f.write("{0} {1} {2}\n".format(dMu2, dMu3, dMu4))
        f.write("{0} {1} {2}\n".format(dMu2Err, dMu3Err, dMu4Err))

# ax.set_xlim([0.0, KGrid[-1]/Para.kF])
# ax.set_xlabel("$q/k_F$", size=size)
# ax.set_ylabel("$-P(\\tau=0^-, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

# plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
# plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
# plt.show()
