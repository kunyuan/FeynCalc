#!/usr/bin/python3
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
Quan = 'LFC'
plt_type = 0  #0: (optimzed) G(q) order 5,  1: G(q) with different lambda

if plt_type ==0:
    fig, ax1 = plt.subplots()
elif plt_type ==1:
    fig, ax = plt.subplots(1,7,sharey=True)
points = ['o','^','s','v','p','<','h','>']

D = 3
Spin = 2

Gfactor = np.zeros(10)
Gfactor_Err = np.zeros(10)
# inlist = open("./inlist", 'r')
inlist = open("./paras", 'r')
for index, eachline in enumerate(inlist):
    para = eachline.split()
    if(len(para)==0):
        break
    beta = float(para[1])
    rs   = float(para[2])
    lam  = float(para[4])

    # with open('../PIMC_data/PIMC_rs{0}_theta{1:3.1f}.txt'.format(rs,1/beta), 'r') as file:
    #     PIMC_para = file.readline()
    #     G_PIMC= np.loadtxt(file)

    with open("./parameter", "w") as file:
        file.write(eachline+"\n")
        file.write("#Order, Beta, rs, Mass2, Lambda, MaxExtMom(*kF), TotalStep(*1e6), Seed, PID")

    Para = param(D, Spin)
    # load chemical potential shift from the file
    mu = np.loadtxt(Para.DataFolder+"/dMu_beta{0}_rs{1}_lam{2}".format(beta,rs,lam))
    dMu2, dMu3, dMu4 = mu[0, :]
    dMu2Err, dMu3Err, dMu4Err = mu[1, :]

    DataDict, Step, Groups, ReWeight, Grids = LoadFile(Para.DataFolder, "pid[0-9]+.dat")
    KGrid = Grids["KGrid"]
    # print (Groups)

    ###### Calculate finite-temperature polarization ################
    BubbleQ = np.zeros((len(KGrid),2))
    for qi, q in enumerate(KGrid):
        BubbleQ[qi] = bubble.Bubble(D, Para.Beta, Spin, Para.kF, q)
        # print ("{0:10.6f}  {1:10.6f}".format(q/Para.kF, BubbleQ[qi]))
        # print (q/Para.kF, BubbleQ[qi][0], BubbleQ[qi][1])
    BubbleQ[0,0] = -BubbleQ[0,0]    

    Bubble = bubble.Bubble(D, Para.Beta, Spin, Para.kF, 0.0)
    print ("Uniform Polarization: ", Bubble[0], "+-", Bubble[1])
    print ("Uniform polarization for the Free electron at T=0: ", Para.Nf)
    Phys = Bubble[0]*len(KGrid)

    EsDataDict = {}
    for g in DataDict.keys():
        # print (g)
        EsDataDict[g] = reduce.EstimateGroup(DataDict, Step, Phys, g)

    for (o, key) in enumerate(sorted(EsDataDict.keys())):
        if key == (0, ):
            continue
        y = EsDataDict[key]
        print (yellow("(Order: {0}, VerCT: {1}, SigmaCT: {2}) = {3:12.8f} +- {4:12.8f}".format(
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
        print (green("(Order: {0}, SigmaCT: {1}) = {2:12.8f} +- {3:12.8f}".format(
            key[0], key[1], y[0][0], y[1][0]*2.0)))

    #print(EsData)

    Each = {}
    if Para.Order >= 1:
        Each[1] = EsData[(1, 0)]
    if Para.Order >= 2:
        Each[2] = EsData[(2, 0)]
    if Para.Order >= 3:
        Each[3] = EsData[(3, 0)]+dMu2*EsData[(1, 1)]
    if Para.Order >= 4:
        Each[4] = EsData[(4, 0)]+dMu2*EsData[(2, 1)]+dMu3*EsData[(1, 1)]
    if Para.Order >= 5:
        Each[5] = EsData[(5, 0)]+dMu2*EsData[(3, 1)]+dMu3 * \
            EsData[(2, 1)]+dMu4*EsData[(1, 1)]+dMu2**2*EsData[(1, 2)]

    Accu = {}
    for o in range(1, Para.Order+1):
        Accu[o] = np.array(Each[o])
        for i in range(1, o):
            Accu[o] += Each[i]

    for o in range(1, Para.Order+1):
        print ("Order {0}: {1:12.8f} +-{2:12.8f}, Accu: {3:12.8f} +-{4:12.8f}".format(
        o, Each[o][0, 0], Each[o][1, 0]*1.0, Accu[o][0, 0], Accu[o][1, 0]*1.0))

    ###### Calculate static local field factor ################
    G = {}
    G_Err = {} 
    for o in range(3, Para.Order+1):
        G[o] = KGrid**2/8.0/np.pi*(-1.0/BubbleQ[:, 0]-1.0/Accu[o][0, :])
        G_Err[o] = KGrid**2/8.0/np.pi* np.sqrt(BubbleQ[:, 1]**2.0/BubbleQ[:, 0]**4.0 +Accu[o][1, :]**2.0/Accu[o][0, :]**4.0)

    G_Err[Para.Order] = G_Err[Para.Order] + np.abs(G[Para.Order]-G[Para.Order-1])

    if plt_type==1:
        ax1 = ax[index]
        omin = 3
    elif plt_type == 0:
        omin = 5
        # G_Err[Para.Order] = G_Err[Para.Order] + np.abs(G[Para.Order] - G[Para.Order-1])
        for qi, q in enumerate(KGrid):
            if qi == 9:
                print(q, Accu[o][0, qi], G[Para.Order][qi])

            if Gfactor[qi]<G[Para.Order][qi]:
                Gfactor[qi] = G[Para.Order][qi]
                Gfactor_Err[qi] = G_Err[Para.Order][qi]

    # print(Gfactor)
    # plt.figure()
    if plt_type == 0:
        ax1.errorbar(KGrid/Para.kF, G[Para.Order], yerr=G_Err[Para.Order], fmt='o-', marker=points[index], capthick=1, capsize=1, ms=6,
                    color=ColorList[index], label=r"Order {0}, $\lambda={1:3.1f}$".format(o,lam))
        print ('Plot G(q)')
    elif plt_type == 1:
        for o in range(omin, Para.Order+1):
            ax1.errorbar(KGrid/Para.kF, G[o], yerr=G_Err[o], fmt='o-', marker=points[0], capthick=1, capsize=1, ms=6,
                        color=ColorList[o-omin+index*3], label=r"Order {0}, $\lambda={1:3.1f}$".format(o,lam))

    if index==0:
        ax1.set_ylabel("$G(\omega=0, q)$", size=size)
        # ax1.set_ylim([0.0, 1.4])
    ax.set_xlim([0.0, KGrid[-1]/Para.kF])
    ax1.set_xlabel("$q/k_F$", size=size)
    ax1.legend(loc=2, frameon=False, fontsize=size)

    # ax1.errorbar(G_PIMC[:,0], G_PIMC[:,1], yerr=G_PIMC[:,2], fmt='o-', marker='3', capthick=1, capsize=1, ms=6,
    #             color=ColorList[-2], label='PIMC')

plt.suptitle(r"$\beta={0}, r_s={1}$".format(beta,rs))
# plt.tight_layout()
if plt_type==0:
    with open('../PIMC_data/PIMC_rs{0}_theta{1:3.1f}.txt'.format(rs,1/beta), 'r') as file:
        PIMC_para = file.readline()
        G_PIMC= np.loadtxt(file)
    ax1.errorbar(G_PIMC[:,0], G_PIMC[:,1], yerr=G_PIMC[:,2], fmt='o-', marker='3', capthick=1, capsize=1, ms=6,
                color=ColorList[-2], label='PIMC')
    ax1.errorbar(KGrid/Para.kF, Gfactor, yerr=Gfactor_Err, fmt='o-', marker='o', capthick=1, capsize=1, ms=6,
                color= ColorList[-1], label="optimzied $G(q)$")
    ax1.legend(loc=2, frameon=False, fontsize=size)
    plt.savefig('Plot/'+Quan+'_charge_beta{0}_rs{1}_opt.pdf'.format(beta,rs))
elif plt_type==1:
    fig.subplots_adjust(wspace=0)
    plt.savefig('Plot/'+Quan+'_charge_beta{0}_rs{1}_lam.pdf'.format(beta,rs))
plt.show()
