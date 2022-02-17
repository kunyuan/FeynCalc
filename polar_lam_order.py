#!/usr/bin/python3
from IO import *
# import scipy.integrate as integrate
# import bubble_dynamic
import bubble
import reduce
# import matplotlib as mat
# import matplotlib.pyplot as plt


D = 3
Spin = 2

list_polar = []
polar_diag = []

inlist = open("./inlist", 'r')
for index, eachline in enumerate(inlist):
    para = eachline.split()
    if(len(para) == 0):
        break
    beta = float(para[1])
    rs = float(para[2])
    lam = float(para[4])

    with open("./parameter", "w") as file:
        file.write(' '.join(para[:-2])+"\n\n")
        file.write(
            "#Order, Beta, rs, Mass2, Lambda, MinExtMom(*kF), MaxExtMom, TotalStep(*1e6)")

    Para = param(D, Spin)
    # load chemical potential shift from the file
    if Para.Order >= 3:
        mu = np.loadtxt(Para.DataFolder +
                        "/dMu_beta{0}_rs{1}_lam{2}".format(beta, rs, lam))
        dMu2, dMu3, dMu4 = mu[0, :]
        dMu2Err, dMu3Err, dMu4Err = mu[1, :]

    DataDict, Step, Groups, ReWeight, Grids = LoadFile(
        Para.DataFolder, "pid[0-9]+.dat")
    KGrid = Grids["KGrid"]
    print(Groups)

    ###### Calculate finite-temperature polarization ################
    BubbleQ = np.zeros(len(KGrid))
    for qi, q in enumerate(KGrid):
        BubbleQ[qi] = bubble.Bubble(D, Para.Beta, Spin, Para.kF, q)[0]
        # print ("{0:10.6f}  {1:10.6f}".format(q/Para.kF, BubbleQ[qi]))
    BubbleQ[0] = -BubbleQ[0]

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
        y = EsDataDict[key]
        print(yellow("(Order: {0}, VerCT: {1}, SigmaCT: {2}) = {3:12.8f} +- {4:12.8f}".format(
            key[0], key[1], key[2], y[0][0], y[1][0]*2.0)))
        EsDataDict[key] = np.array((y[0], y[1]))
        for qi, q in enumerate(KGrid):
            dat = np.array(
                [q, y[0][qi], y[1][qi], int(key[0]*100+key[1]*10+key[2])])
            polar_diag.append(dat)

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

    # print(EsData)

    Each = {}
    if Para.Order >= 1:
        Each[1] = EsData[(1, 0)]
    if Para.Order >= 2:
        Each[2] = EsData[(2, 0)]
    if Para.Order >= 3:
        Each[3] = EsData[(3, 0)]+dMu2*EsData[(1, 1)]
        # Each[3][1,:] = dMu2Err**2.0*EsData[(1, 1)][0,:]**2.0 + dMu2**2.0*EsData[(1, 1)][1,:]**2.0
        Each[3][1, :] = dMu2**2.0*EsData[(1, 1)][1, :]**2.0
        Each[3][1, :] = (Each[3][1, :] + EsData[(3, 0)][1, :]**2.0)**0.5

    if Para.Order >= 4:
        Each[4] = EsData[(4, 0)]+dMu2*EsData[(2, 1)]+dMu3*EsData[(1, 1)]
        # Each[4][1,:] = dMu2Err**2.0*EsData[(2, 1)][0,:]**2.0 + dMu2**2.0*EsData[(2, 1)][1,:]**2.0
        # Each[4][1,:]+= dMu3Err**2.0*EsData[(1, 1)][0,:]**2.0 + dMu3**2.0*EsData[(1, 1)][1,:]**2.0
        Each[4][1, :] = dMu2**2.0 * \
            EsData[(2, 1)][1, :]**2.0 + dMu3**2.0*EsData[(1, 1)][1, :]**2.0
        Each[4][1, :] = (Each[4][1, :] + EsData[(4, 0)][1, :]**2.0)**0.5
    if Para.Order >= 5:
        Each[5] = EsData[(5, 0)]+dMu2*EsData[(3, 1)]+dMu3 * \
            EsData[(2, 1)]+dMu4*EsData[(1, 1)]+dMu2**2*EsData[(1, 2)]
        # Each[5][1,:] = dMu2Err**2.0*EsData[(3, 1)][0,:]**2.0 + dMu2**2.0*EsData[(3, 1)][1,:]**2.0
        # Each[5][1,:]+= dMu3Err**2.0*EsData[(2, 1)][0,:]**2.0 + dMu3**2.0*EsData[(2, 1)][1,:]**2.0
        # Each[5][1,:]+= dMu4Err**2.0*EsData[(1, 1)][0,:]**2.0 + dMu4**2.0*EsData[(1, 1)][1,:]**2.0
        # Each[5][1,:]+= 4*dMu2Err**2.0*(dMu2*EsData[(1, 2)][0,:])**2.0 + dMu2**4.0*EsData[(1, 2)][1,:]**2.0
        Each[5][1, :] = dMu2**2.0*EsData[(3, 1)][1, :]**2.0 + dMu3**2.0*EsData[(2, 1)][1, :]**2.0 + \
            dMu4**2.0*EsData[(1, 1)][1, :]**2.0 + dMu2**4.0 * \
            EsData[(1, 2)][1, :]**2.0
        Each[5][1, :] = (Each[5][1, :] + EsData[(5, 0)][1, :]**2.0)**0.5

    Accu = {}
    for o in range(1, Para.Order+1):
        Accu[o] = np.array(Each[o])
        for i in range(1, o):
            Accu[o] += Each[i]
    # print Each[1].shape

    for o in range(1, Para.Order+1):
        print("Order {0}: {1:12.8f} +-{2:12.8f}, Accu: {3:12.8f} +-{4:12.8f}".format(
            o, Each[o][0, 0], Each[o][1, 0]*1.0, Accu[o][0, 0], Accu[o][1, 0]*1.0))
        for qi, q in enumerate(KGrid):
            # dat = np.array([q, Accu[o][0, qi], Accu[o][1, qi], o, lam, beta, rs])
            dat = np.array([q, Each[o][0, qi], Each[o]
                           [1, qi], o, lam, beta, rs])
            list_polar.append(dat)

BubbleQ = np.zeros((len(KGrid), 2))
polar0 = []
for qi, q in enumerate(KGrid):
    BubbleQ[qi] = bubble.Bubble(D, Para.Beta, Spin, Para.kF, q)
    if qi == 0:
        BubbleQ[0][0] = -BubbleQ[0][0]
    dat = np.array([q, BubbleQ[qi, 0], BubbleQ[qi, 1]])
    polar0.append(dat)
    # print ("{0:10.6f}  {1:10.6f}".format(q/Para.kF, BubbleQ[qi]))
polar0 = np.array(polar0)
with open('./polar0_beta{0:3.1f}_rs{1:3.1f}.txt'.format(beta, rs), 'w') as file:
    file.write("# q Polar Error\n")
    np.savetxt(file, polar0, delimiter=" ")

list_polar = np.array(list_polar)
with open('./Data_Polar/polarO_beta{0:3.1f}_rs{1:3.1f}_o{2}.dat'.format(beta, rs, Para.Order), 'w') as file:
    file.write("# q Polar Error order lambda beta rs\n")
    np.savetxt(file, list_polar, delimiter=" ")

polar_diag = np.array(polar_diag)

with open('./Data_Polar/polarDiag_beta{0:3.1f}_rs{1:3.1f}_o{2}.dat'.format(beta, rs, Para.Order), 'w') as file:
    file.write("# q Polar Error DiagID\n")
    np.savetxt(file, polar_diag, delimiter=" ")

    # ###### Calculate static local field factor ################
    # G = {}
    # for o in range(1, Para.Order+1):
    #     G[o] = KGrid**2/8.0/np.pi*(-1.0/BubbleQ-1.0/Accu[o][0, :])
    #     print(Accu[o][0,13], Accu[o][1,13]*2.0)
