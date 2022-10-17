#!/usr/bin/python3
from IO import *
# import bubble_dynamic
import bubble
import reduce

D = 3
Spin = 2

polarall = []
polar0 = []
inlist = open("./paras", 'r')
for index, eachline in enumerate(inlist):
    para = eachline.split()
    if(len(para) == 0):
        break
    beta = float(para[1])
    rs = float(para[2])
    lam = float(para[4])

    with open("./parameter", "w") as file:
        file.write(eachline+"\n")
        file.write(
            "#Order, Beta, rs, Mass2, Lambda, MaxExtMom(*kF), TotalStep(*1e6), Seed, PID")

    Para = param(D, Spin)
    # load chemical potential shift from the file
    if Para.Order >= 3:
        mu = np.loadtxt(
            "../dMu/dMu_beta{0}_rs{1}_lam{2}".format(beta, rs, lam))
        dMu2, dMu3, dMu4 = mu[0, :]
        dMu2Err, dMu3Err, dMu4Err = mu[1, :]

    DataDict, Step, Groups, ReWeight, Grids = LoadFile(
        Para.DataFolder, "pid[0-9]+.dat")
    KGrid = Grids["KGrid"]
    # print (Groups)
    Iflq = Para.DataFolder[-2:] == 'lq'
    Ifhq = Para.DataFolder[-2:] == 'hq'

    ###### Calculate finite-temperature polarization ################
    BubbleQ = np.zeros((len(KGrid), 2))
    BubbleQCT = np.zeros((len(KGrid), 2))
    for qi, q in enumerate(KGrid):
        BubbleQ[qi] = bubble.Bubble(D, Para.Beta, Spin, Para.kF, q)
        BubbleQCT[qi] = bubble.BubbleCT(D, Para.Beta, Spin, Para.kF, q)
        print("{0:3.1f}".format(q/Para.kF), BubbleQCT[qi][0], BubbleQCT[qi][1])
        # print ("{0:10.6f}  {1:10.6f}".format(q/Para.kF, BubbleQ[qi]))
    # BubbleQ[0] = -BubbleQ[0]

    # print (BubbleQ)

    Bubble = bubble.Bubble(D, Para.Beta, Spin, Para.kF, 0.0)
    print("Uniform Polarization: ", Bubble[0], "+-", Bubble[1])
    print("Uniform polarization for the Free electron at T=0: ", Para.Nf)
    Phys = Bubble[0]*len(KGrid)

    EsDataDict = {}
    for g in DataDict.keys():
        print(g)
        EsDataDict[g] = reduce.EstimateGroup(DataDict, Step, Phys, g)

    for (o, key) in enumerate(sorted(EsDataDict.keys())):
        if key == (0, ):
            continue
        # if o == 1:
        #     EsDataDict[key][0] += BubbleQ[:, 0]
        # if key == (1, 0, 1, 2):
        #     EsDataDict[key][0] += BubbleQCT[:, 0]
        y = EsDataDict[key]
        print(yellow("(Order: {0}, VerCT: {1}, SigmaCT: {2}) = {3:12.8f} +- {4:12.8f}".format(
            key[0], key[1], key[2], y[0][0], y[1][0]*2.0)))
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
        print(green("(Order: {0}, SigmaCT: {1}) = {2:12.8f} +- {3:12.8f}".format(
            key[0], key[1], y[0][0], y[1][0]*1.0)))

    Each = {}
    if Para.Order >= 1:
        Each[1] = EsData[(1, 0)]
    if Para.Order >= 2:
        Each[2] = EsData[(2, 0)]
    if Para.Order >= 3:
        Each[3] = EsData[(3, 0)]+dMu2*EsData[(1, 1)]
        Each[3][1, :] = dMu2Err**2.0 * \
            EsData[(1, 1)][0, :]**2.0 + dMu2**2.0*EsData[(1, 1)][1, :]**2.0
        Each[3][1, :] = (Each[3][1, :] + EsData[(3, 0)][1, :]**2.0)**0.5
        # print(dMu2Err*EsData[(1, 1)][0,:])
        # print(dMu2*EsData[(1, 1)][1,:])
        # print(EsData[(3, 0)][1,:])

    if Para.Order >= 4:
        Each[4] = EsData[(4, 0)]+dMu2*EsData[(2, 1)]+dMu3*EsData[(1, 1)]
        Each[4][1, :] = dMu2Err**2.0 * \
            EsData[(2, 1)][0, :]**2.0 + dMu2**2.0*EsData[(2, 1)][1, :]**2.0
        Each[4][1, :] += dMu3Err**2.0 * \
            EsData[(1, 1)][0, :]**2.0 + dMu3**2.0*EsData[(1, 1)][1, :]**2.0
        Each[4][1, :] = (Each[4][1, :] + EsData[(4, 0)][1, :]**2.0)**0.5
    if Para.Order >= 5:
        Each[5] = EsData[(5, 0)]+dMu2*EsData[(3, 1)]+dMu3 * \
            EsData[(2, 1)]+dMu4*EsData[(1, 1)]+dMu2**2*EsData[(1, 2)]
        Each[5][1, :] = dMu2Err**2.0 * \
            EsData[(3, 1)][0, :]**2.0 + dMu2**2.0*EsData[(3, 1)][1, :]**2.0
        Each[5][1, :] += dMu3Err**2.0 * \
            EsData[(2, 1)][0, :]**2.0 + dMu3**2.0*EsData[(2, 1)][1, :]**2.0
        Each[5][1, :] += dMu4Err**2.0 * \
            EsData[(1, 1)][0, :]**2.0 + dMu4**2.0*EsData[(1, 1)][1, :]**2.0
        Each[5][1, :] += 4*dMu2Err**2.0 * \
            (dMu2*EsData[(1, 2)][0, :])**2.0 + \
            dMu2**4.0*EsData[(1, 2)][1, :]**2.0
        Each[5][1, :] = (Each[5][1, :] + EsData[(5, 0)][1, :]**2.0)**0.5

    Accu = {}
    for o in range(1, Para.Order+1):
        Accu[o] = np.array(Each[o])
        for i in range(1, o):
            Accu[o] += Each[i]

    list_polar = []
    for o in range(1, Para.Order+1):
        # print("Order {0}: {1:12.8f} +-{2:12.8f}, Accu: {3:12.8f} +-{4:12.8f}".format(
        #     o, Each[o][0, 0], Each[o][1, 0]*1.0, Accu[o][0, 0], Accu[o][1, 0]*1.0))
        for qi, q in enumerate(KGrid):
            dat = np.array([q, Accu[o][0, qi], Accu[o]
                           [1, qi], o, lam, beta, rs])
            list_polar.append(dat)
            polarall.append(dat)

    BubbleQ = np.zeros((len(KGrid), 2))

    for qi, q in enumerate(KGrid):
        BubbleQ[qi] = bubble.Bubble(D, Para.Beta, Spin, Para.kF, q)
        if qi == 0:
            BubbleQ[0][0] = -BubbleQ[0][0]
        dat = np.array([q, BubbleQ[qi, 0], BubbleQ[qi, 1]])
        polar0.append(dat)
        print("q = {0}kF".format(q/Para.kF))
        for o in range(1, Para.Order+1):
            print("Order {0}: {1:12.8f} +-{2:12.8f}, Accu: {3:12.8f} +-{4:12.8f}".format(
                o, Each[o][0, qi], Each[o][1, qi]*1.0, Accu[o][0, qi], Accu[o][1, qi]*1.0))

    list_polar = np.array(list_polar)
    filename = './Data_Polar/polar_beta{0}_rs{1:3.1f}_lam{2:.1f}_o{3}.dat'.format(
        beta, rs, lam, Para.Order)
    with open(filename, 'a') as file:
        print('Save to '+filename)
        file.write("# q Polar Error order lambda beta rs\n")
        np.savetxt(file, list_polar, delimiter=" ")

polarall = np.array(polarall)
# filename = './Data_Polar/polar_beta{0}_rs{1:3.1f}_o{2}.dat'.format(beta, rs, Para.Order)
filename = './Data_Polar/polar_beta{0}_rs{1:3.1f}_o{2}.dat'.format(
    beta, rs, Para.Order)
with open(filename, 'a') as file:
    print('Save to '+filename)
    file.write("# q Polar Error order lambda beta rs\n")
    np.savetxt(file, polarall, delimiter=" ")

polar0 = np.array(polar0)
polar0 = polar0[polar0[:, 0].argsort()]
with open('./Data_Polar/polar0_beta{0}_rs{1:3.1f}.txt'.format(beta, rs), 'a') as file:
    # file.write("# q Polar Error\n")
    np.savetxt(file, polar0, delimiter=" ")
