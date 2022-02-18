# Calculate kinetic energy
# !/usr/bin/python
from IO import *
import bubble
import kinetic
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
Kinetic = []
for index, eachline in enumerate(inlist):
    para = eachline.split()
    if len(para) == 0:
        break
    beta = float(para[1])
    rs = float(para[2])
    lam = float(para[4])
    with open("./parameter", "w") as file:
        parameters = ' '.join(para[:-2])
        file.write(parameters+"\n\n")
        file.write(
            "#Order, Beta, rs, Mass2, Lambda, MinExtMom(*kF), MaxExtMom, TotalStep(*1e6)")

    Para = param(D, Spin)

    mu = np.loadtxt("./dMu/dMu_beta{0}_rs{1}_lam{2}".format(beta, rs, lam))
    dMu2, dMu3, dMu4 = mu[0, :]
    dMu2Err, dMu3Err, dMu4Err = mu[1, :]

    DataDict, Step, Groups, ReWeight, Grids = LoadFile(
        Para.DataFolder, "Ek_pid[0-9]+.dat")
    KGrid = Grids["KGrid"]
    # for i in range(len(Data)):
    #     print i, Data[i][1, 0]

    Bubble = bubble.Bubble(D, Para.Beta, Spin, Para.kF, 0.0)
    Ek0 = kinetic.Kinetic(D, Para.Beta, Spin, Para.kF)
    print("Noninteracting kinetic energy: ", Ek0[0], "+-", Ek0[1])
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
        if o == 1:
            EsDataDict[key][0] += Ek0[0]
        y = EsDataDict[key] / Para.Beta
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

    Each = {}
    if Para.Order >= 1:
        Each[1] = EsData[(1, 0)]
    if Para.Order >= 2:
        # Each[2] = EsData[(2, 0)]
        Each[2] = [0.0, 0.0]
    if Para.Order >= 3:
        Each[3] = EsData[(3, 0)]+dMu2*EsData[(1, 1)]
        Each[3][1] = dMu2Err**2.0 * \
            EsData[(1, 1)][0]**2.0 + dMu2**2.0*EsData[(1, 1)][1]**2.0
        Each[3][1] = dMu2**2.0*EsData[(1, 1)][1]**2.0
        Each[3][1] = (Each[3][1] + EsData[(3, 0)][1]**2.0)**0.5
    if Para.Order >= 4:
        Each[4] = EsData[(4, 0)]+dMu2*EsData[(2, 1)]+dMu3*EsData[(1, 1)]
        Each[4][1] = dMu2Err**2.0 * \
            EsData[(2, 1)][0]**2.0 + dMu2**2.0*EsData[(2, 1)][1]**2.0
        Each[4][1] += dMu3Err**2.0 * \
            EsData[(1, 1)][0]**2.0 + dMu3**2.0*EsData[(1, 1)][1]**2.0
        Each[4][1] = (Each[4][1] + EsData[(4, 0)][1]**2.0)**0.5
    if Para.Order >= 5:
        Each[5] = EsData[(5, 0)]+dMu2*EsData[(3, 1)]+dMu3 * \
            EsData[(2, 1)]+dMu4*EsData[(1, 1)]+dMu2**2*EsData[(1, 2)]
        Each[5][1] = dMu2Err**2.0 * \
            EsData[(3, 1)][0]**2.0 + dMu2**2.0*EsData[(3, 1)][1]**2.0
        Each[5][1] += dMu3Err**2.0 * \
            EsData[(2, 1)][0]**2.0 + dMu3**2.0*EsData[(2, 1)][1]**2.0
        Each[5][1] += dMu4Err**2.0 * \
            EsData[(1, 1)][0]**2.0 + dMu4**2.0*EsData[(1, 1)][1]**2.0
        Each[5][1] += 4*dMu2Err**2.0 * \
            (dMu2*EsData[(1, 2)][0])**2.0 + dMu2**4.0*EsData[(1, 2)][1]**2.0
        Each[5][1] = (Each[5][1] + EsData[(5, 0)][1]**2.0)**0.5
    # plt.errorbar(KGrid/Para.kF, y, yerr=err, fmt='o-', capthick=1, capsize=4,
    #              color=ColorList[o], label="Order: {0}, SigamCT: {1}".format(key[0], key[1]))
    Accu = {}
    for o in range(1, Para.Order+1):
        Accu[o] = np.array(Each[o])
        Accu[o][1] = Each[o][1]**2.0
        for i in range(1, o):
            Accu[o][0] += Each[i][0]
            Accu[o][1] += Each[i][1]**2.0
        Accu[o][1] = Accu[o][1]**0.5
        print("Order {0}: {1:12.8f} +-{2:12.8f}, Kinetic_Accu: {3:12.8f} +-{4:12.8f}".format(
            o, Each[o][0], Each[o][1]*1.0, Accu[o][0], Accu[o][1]))

    density = 3.0/(4*np.pi * rs**3)
    print('Density * temperature = ', density*Para.Beta)
    filename = "./kinetic.dat"
    with open(filename, "a") as f:
        print('Save Kinetic energy data in '+filename)
        f.write("{0} {1} {2} {3} {4}\n".format(
            Accu[Para.Order][0], Accu[Para.Order][1], beta, rs, lam))
