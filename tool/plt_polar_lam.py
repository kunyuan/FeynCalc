#!/usr/bin/python3
from IO import *
import scipy.integrate as integrate
#import bubble_dynamic
import bubble
import reduce
import matplotlib as mat
import matplotlib.pyplot as plt

plt.switch_backend('TkAgg')
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12
Quan = 'LFC'
plt_type = 0     #0: (optimzed) G(q) order 5,  1: G(q) with fixQ
fig, ax1 = plt.subplots()

points = ['o','^','s','v','p','<','h','>']
pltindex = {0.125:0, 0.25:1, 0.5:2, 1.0:3, 2.0:4, 4.0:5, 8.0:6, 16.0:7}

D = 3
Spin = 2

kmax_G = {0.125: 6.2, 0.25: 4.4, 0.5: 3.1, 1.0: 2.5, 2.0:2.3, 4.0:2, 8.0:2, 16.0:2}   #<1: 2.2/sqrt(beta); >2: 2
Gfactor = {}
Gfixq = {}
Piq4_all = {}
Piq4fixq = {}
fixQ = [2.0, 3.0, 4.0, 4.5, 5.0]
inlist = open("./paras", 'r')
for index, eachline in enumerate(inlist):
    para = eachline.split()
    if(len(para)==0):
        break
    beta = float(para[1])
    rs   = float(para[2])
    lam  = float(para[4])
    qmax = float(para[5])

    # with open('../PIMC_data/PIMC_rs{0}_theta{1:3.1f}.txt'.format(rs,1/beta), 'r') as file:
    #     PIMC_para = file.readline()
    #     G_PIMC= np.loadtxt(file)

    with open("./parameter", "w") as file:
        file.write(eachline+"\n")
        file.write("#Order, Beta, rs, Mass2, Lambda, MaxExtMom(*kF), TotalStep(*1e6), Seed, PID")

    Para = param(D, Spin)
    # load chemical potential shift from the file
    if Para.Order >= 3:
        mu = np.loadtxt(Para.DataFolder+"/dMu_beta{0}_rs{1}_lam{2}".format(beta,rs,lam))
        dMu2, dMu3, dMu4 = mu[0, :]
        dMu2Err, dMu3Err, dMu4Err = mu[1, :]

    # load polarization of order 1
    fname = "./Data_Polar/polar_beta{0}_rs{1}_o1.dat".format(beta,rs)
    IFinp1 = False
    try:
        if qmax == 10.0:
            with open(fname, "r") as file:
                print ("Loading ", fname)
                polar1 = np.loadtxt(file)
                IFinp1 = True
    except:
        IFinp1 = False
    # IFinp1 = False

    # DataDict, Step, Groups, ReWeight, Grids = LoadFile(Para.DataFolder, "pid[0-9]+.dat")
    DataDict, Step, Groups, ReWeight, Grids = LoadFile(Para.DataFolder, "pid")
    KGrid = Grids["KGrid"]
    # print (Groups)

    ###### Calculate finite-temperature polarization ################
    BubbleQ = np.zeros((len(KGrid),2))
    for qi, q in enumerate(KGrid):
        BubbleQ[qi] = bubble.Bubble(D, Para.Beta, Spin, Para.kF, q)  
        # print ("{0:3.1f}".format(q/Para.kF), BubbleQ[qi][0], BubbleQ[qi][1])
    # BubbleQ[0,0] = -BubbleQ[0,0]    

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
            key[0], key[1], key[2], y[0][0], y[1][0])))
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
            key[0], key[1], y[0][0], y[1][0])))

    # if beta == 1.0:
    #     print('beta 1  !')
    #     Para.Order = 3

    Each = {}
    if Para.Order >= 1:
        # print(EsData[(1, 0)])
        if IFinp1:
            # print(polar1[:,1],polar1[:,2])
            Each[1] = np.array([polar1[:,1],polar1[:,2]])
            print (red("(Order: 1, SigmaCT: 0) = {0:12.8f} +- {1:12.8f}".format(polar1[0][1], polar1[0][2])))
        else:
            Each[1] = EsData[(1, 0)]
    if Para.Order >= 2:
        Each[2] = EsData[(2, 0)]
        # print(Each[2])
    if Para.Order >= 3:
        Each[3] = EsData[(3, 0)]+dMu2*EsData[(1, 1)]
        Each[3][1,:] = dMu2Err**2.0*EsData[(1, 1)][0,:]**2.0 + dMu2**2.0*EsData[(1, 1)][1,:]**2.0 
        Each[3][1,:] = dMu2**2.0*EsData[(1, 1)][1,:]**2.0 
        Each[3][1,:] = (Each[3][1,:] + EsData[(3, 0)][1,:]**2.0)**0.5
        print("order 3")
        print(EsData[(3, 0)][1,:])         ## No.2
        print(dMu2Err*EsData[(1, 1)][0,:]) ## No.1
        print(dMu2*EsData[(1, 1)][1,:])
        # print(Each[3])
    if Para.Order >= 4:
        Each[4] = EsData[(4, 0)]+dMu2*EsData[(2, 1)]+dMu3*EsData[(1, 1)]
        Each[4][1,:] = dMu2Err**2.0*EsData[(2, 1)][0,:]**2.0 + dMu2**2.0*EsData[(2, 1)][1,:]**2.0 
        Each[4][1,:]+= dMu3Err**2.0*EsData[(1, 1)][0,:]**2.0 + dMu3**2.0*EsData[(1, 1)][1,:]**2.0 
        Each[4][1,:] = (Each[4][1,:] + EsData[(4, 0)][1,:]**2.0)**0.5
        print("order 4")
        print(EsData[(4, 0)][1,:])         ## No.2
        print(dMu2Err*EsData[(2, 1)][0,:])
        print(dMu2*EsData[(2, 1)][1,:])
        print(dMu3Err*EsData[(1, 1)][0,:]) ## No.1
        print(dMu3*EsData[(1, 1)][1,:])
    if Para.Order >= 5:
        Each[5] = EsData[(5, 0)]+dMu2*EsData[(3, 1)]+dMu3 * \
            EsData[(2, 1)]+dMu4*EsData[(1, 1)]+dMu2**2*EsData[(1, 2)]
        Each[5][1,:] = dMu2Err**2.0*EsData[(3, 1)][0,:]**2.0 + dMu2**2.0*EsData[(3, 1)][1,:]**2.0 
        Each[5][1,:]+= dMu3Err**2.0*EsData[(2, 1)][0,:]**2.0 + dMu3**2.0*EsData[(2, 1)][1,:]**2.0 
        Each[5][1,:]+= dMu4Err**2.0*EsData[(1, 1)][0,:]**2.0 + dMu4**2.0*EsData[(1, 1)][1,:]**2.0 
        Each[5][1,:]+= 4*dMu2Err**2.0*(dMu2*EsData[(1, 2)][0,:])**2.0 + dMu2**4.0*EsData[(1, 2)][1,:]**2.0 
        Each[5][1,:] = (Each[5][1,:] + EsData[(5, 0)][1,:]**2.0)**0.5
        print("order 5")
        print(EsData[(5, 0)][1,:])         ## No.2
        print(dMu2Err*EsData[(3, 1)][0,:])
        print(dMu2*EsData[(3, 1)][1,:])
        print(dMu3Err*EsData[(2, 1)][0,:])
        print(dMu3*EsData[(2, 1)][1,:])
        print(dMu4Err*EsData[(1, 1)][0,:]) ## No.1
        print(dMu4*EsData[(1, 1)][1,:])
        print(2* dMu2Err* dMu2* EsData[(1, 2)][0,:])
        print(dMu2**2 *EsData[(1, 2)][1,:])
    Accu = {}
    for o in range(1, Para.Order+1):
        Accu[o] = np.array(Each[o])
        Accu[o][1,:] = Each[o][1,:]**2.0
        for i in range(1, o):
            Accu[o][0,:] += Each[i][0,:]
            Accu[o][1,:] += Each[i][1,:]**2.0
        Accu[o][1,:] = Accu[o][1,:]**0.5

    for o in range(1, Para.Order+1):
        print ("Order {0}: {1:12.8f} +-{2:12.8f}, Accu: {3:12.8f} +-{4:12.8f}".format(
        o, Each[o][0, 0], Each[o][1, 0]*1.0, Accu[o][0, 0], Accu[o][1, 0]*1.0))

    ###### Calculate static local field factor ################
    density = 3.0/2/np.pi/Para.Rs**3.0
    omin = 1
    G = {}
    G_Err = {} 
    Piq4 = {}
    Piq4_Err = {}
    for o in range(omin, Para.Order+1):
        G[o] = KGrid**2/8.0/np.pi*(-1.0/BubbleQ[:, 0]-1.0/Accu[o][0, :])
        G_Err[o] = KGrid**2/8.0/np.pi* np.sqrt(BubbleQ[:, 1]**2.0/BubbleQ[:, 0]**4.0 +Accu[o][1, :]**2.0/Accu[o][0, :]**4.0)
        # G[o] = 1.0/8.0/np.pi*(-1.0/BubbleQ[:, 0]-1.0/Accu[o][0, :])
        # G_Err[o] = 1.0/8.0/np.pi* np.sqrt(BubbleQ[:, 1]**2.0/BubbleQ[:, 0]**4.0 +Accu[o][1, :]**2.0/Accu[o][0, :]**4.0)
        Piq4[o] = KGrid**4/8.0/np.pi*(BubbleQ[:, 0]+Accu[o][0, :])/density**2
        Piq4_Err[o] = KGrid**4/8.0/np.pi* np.abs(Accu[o][1, :])/density**2

    # G_sys = np.abs(G[Para.Order]-G[Para.Order-1])
    # G_Err[Para.Order] = G_Err[Para.Order] + G_sys
    for qi, q in enumerate(KGrid):
        print (q/Para.kF, Accu[o][0 ,qi], Accu[o][1 ,qi], BubbleQ[qi,0], G[Para.Order][qi])
        for qfix in fixQ:
            if abs(q/Para.kF-qfix)<1e-4:
                if qfix in Gfixq.keys():
                    Gfixq[qfix].append([beta, G[Para.Order][qi], G_Err[Para.Order][qi]]) 
                else:
                    Gfixq[qfix] = [[beta, G[Para.Order][qi], G_Err[Para.Order][qi]]]

    # print(G[Para.Order], G_Err[Para.Order])
    # G_lq = G[Para.Order][:,np.newaxis]
    G_lq = np.concatenate((KGrid[:,np.newaxis]/Para.kF,G[Para.Order][:,np.newaxis], G_Err[Para.Order][:,np.newaxis]), axis=1)
    Piq4_lq = np.concatenate((KGrid[:,np.newaxis]/Para.kF,Piq4[Para.Order][:,np.newaxis], Piq4_Err[Para.Order][:,np.newaxis]), axis=1)
    if beta in Gfactor.keys():
        Gfactor[beta] = np.concatenate((Gfactor[beta], G_lq), axis=0)
        Piq4_all[beta] = np.concatenate((Piq4_all[beta], Piq4_lq), axis=0)
    else:
        Gfactor[beta] = G_lq
        Piq4_all[beta] = Piq4_lq
   
    # elif plt_type == 0:
        # G_Err[Para.Order] = G_Err[Para.Order] + np.abs(G[Para.Order] - G[Para.Order-1])
        # SysErr = G_sys
        # for qi, q in enumerate(KGrid):
        #     # if G_sys[qi] <= SysErr[qi]:
        #     #     SysErr[qi] = G_sys[qi]
        #     #     Gfactor[qi] = G[Para.Order][qi]
        #     #     Gfactor_Err[qi] = G_Err[Para.Order][qi]
        #     if Gfactor[qi]<G[Para.Order][qi]:
        #         Gfactor[qi] = G[Para.Order][qi]
        #         Gfactor_Err[qi] = G_Err[Para.Order][qi]

    # if plt_type == 0:
    #     # ax1.errorbar(KGrid/Para.kF, G[Para.Order], yerr=G_Err[Para.Order], fmt='o-', marker=points[pltindex[beta]], capthick=1, capsize=1, ms=6,
    #                 # color=ColorList[pltindex[beta]], label=r"$\beta={0:4.2f},\lambda={1:4.2f}$".format(beta,lam))
    #     print ('Plot G(q)')
    # elif plt_type == 1:
    #     for o in range(omin, Para.Order+1):
    #         ax1.errorbar(KGrid/Para.kF, G[o], yerr=G_Err[o], fmt='o-', marker=points[0], capthick=1, capsize=1, ms=6,
    #                     color=ColorList[o-omin+index*3], label=r"Order {0}, $\lambda={1:4.2f}$".format(o,lam))

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)
    if index==0:
        ax1.set_ylabel("$G(\omega=0, q)$", size=size)
        # ax1.set_ylim([0,1.6])
    # ax1.set_xlim([KGrid[0]/Para.kF-0.1, KGrid[-1]/Para.kF+0.1])

    # ax1.errorbar(G_PIMC[:,0], G_PIMC[:,1], yerr=G_PIMC[:,2], fmt='o-', marker='3', capthick=1, capsize=1, ms=6,
    #             color=ColorList[-2], label='PIMC')

# plt.suptitle(r"$\beta={0}, r_s={1}$".format(beta,rs))
# plt.tight_layout()
if plt_type==0:
    for (idx, key) in enumerate(sorted(Gfactor.keys())):
        ax1.errorbar(Gfactor[key][:,0], Gfactor[key][:,1], yerr=Gfactor[key][:,2], fmt='o-', marker=points[idx], capthick=1, capsize=1, ms=6,
            color= ColorList[idx*2], label=r"$\beta={0:4.3f}$".format(key))
        # ax1.errorbar(Piq4_all[key][:,0], Piq4_all[key][:,1], yerr=Piq4_all[key][:,2], fmt='o-', marker=points[idx], capthick=1, capsize=1, ms=6,
        #     color= ColorList[idx*2+1], label=r"$\beta={0:4.3f}, (\Pi-\Pi_0)q^4/n^2$".format(key))
        # ax1.errorbar(Gfactor[key][:,0]/kmax_G[key], Gfactor[key][:,1], yerr=Gfactor[key][:,2], fmt='o-', marker=points[idx], capthick=1, capsize=1, ms=6,
            # color= ColorList[idx], label=r"$\beta={0:4.3f}$".format(key))
    plt.title(r"Order ${0}, r_s={1}$".format(Para.Order, rs))

    #### new PIMC raw data 
    # PIMCfile = '../PIMC_data/PIMC_rs{0}_theta0.0.txt'.format(rs)
    PIMCfile = '../PIMC_data/PIMC_rs{0}_theta{1:3.1f}.txt'.format(rs,1/beta)
    try:
        with open(PIMCfile, 'r') as file:
            print('Loading '+PIMCfile)
            PIMC_para = file.readline()
            G_PIMC= np.loadtxt(file)
        # ax1.errorbar(G_PIMC[:,0], G_PIMC[:,1], yerr=G_PIMC[:,2], fmt='o-', marker='3', capthick=1, capsize=1, ms=6,
                # color=ColorList[-2], label= 'PIMC') #'GS parametrization') # 
    except:
        pass

    # ### DMC data  at zero temperature
    # fname = '../PIMC_data/DMC_rs{0}_theta0.0.txt'.format(rs)
    # with open(fname, 'r') as file:
    #     print('Loading '+PIMCfile)
    #     G_DMC= np.loadtxt(file)
    # print(G_DMC)
    # G_DMC[:,1] = G_DMC[:,1] *(G_DMC[:,0])**2.0
    # G_DMC[:,2] = G_DMC[:,2] *(G_DMC[:,0])**2.0
    # ax1.errorbar(G_DMC[:,0], G_DMC[:,3], yerr=G_DMC[:,4], fmt='o-', marker='4', capthick=1, capsize=1, ms=6,
    #         color=ColorList[-3], label='DMC')

    ### dielectric data from papar 1994
    # PIMCfile = '../PIMC_data/diele_rs{0}_theta0.0.txt'.format(rs)
    # with open(PIMCfile, 'r') as file:
    #     print('Loading '+PIMCfile)
    #     # PIMC_para = file.readline()
    #     diele= np.loadtxt(file)
    # print(diele)  
    # q_d = diele[:,0]
    # BubbleQ_d = np.zeros((len(q_d),2))
    # for qi, q in enumerate(q_d):
    #     BubbleQ_d[qi] = bubble_static.Bubble(D, 100, Spin, Para.kF, q)
    # G_d = 1.0 + 1.0/(diele[:,3]-1.0) - q_d**2.0/8.0/np.pi/BubbleQ_d[:,0]
    # Err_G_d = np.sqrt(diele[:,4]**2.0/(diele[:,3]-1.0)**4.0 + (q_d**2.0/8.0/np.pi*BubbleQ_d[:,1]/BubbleQ_d[:,0]**2.0)**2.0)
    # ax1.errorbar(q_d/Para.kF, G_d, yerr=Err_G_d, fmt='o-', marker=points[-1], capthick=1, capsize=1, ms=6,
    #             color=ColorList[-4], label=r'from $\epsilon_D$')
    # ax1.set_xlim([0, 3])
    # ax1.set_ylim([-0.02,0.04])
    ax1.set_xlabel("$q/k_F^*$", size=size)
    ax1.legend(loc=3, frameon=False, fontsize=size)
    # plt.savefig('Plot/'+Quan+'_charge_beta{0}_rs{1}_opt.pdf'.format(beta,rs))
    # plt.savefig('Plot/'+Quan+'_charge_rs{0}.pdf'.format(rs))
elif plt_type==1:
    print(Gfixq)
    for (idx, qfix) in enumerate(fixQ):
        # print(Gfixq[qfix])
        Gfixq[qfix] = np.array(Gfixq[qfix])
        print(qfix, Gfixq[qfix])
        ax1.errorbar(Gfixq[qfix][:,0], Gfixq[qfix][:,1], yerr=Gfixq[qfix][:,2], fmt='o-', marker=points[idx], capthick=1, capsize=1, ms=6,
            color= ColorList[idx], label=r"$q={0:3.1f}k_F$".format(qfix))
    # fig.subplots_adjust(wspace=0)
    ax1.set_xlabel(r"$\beta$", size=size)
    plt.xscale('log')
    ax1.legend(loc=4, frameon=False, fontsize=size)
    plt.savefig('Plot/'+Quan+'_charge_rs{0}_fixq{1}.pdf'.format(rs, fixQ))
plt.show()
