#!/usr/bin/python3
from IO import *
# import scipy.integrate as integrate
# import bubble_dynamic
import bubble
import reduce
import matplotlib as mat
import matplotlib.pyplot as plt
from brokenaxes import brokenaxes
from matplotlib.gridspec import GridSpec

plt.switch_backend('TkAgg')
mat.rcParams.update({'font.size': 10})
mat.rcParams["font.family"] = "Times New Roman"
size = 11

D = 3
Spin = 2

ExtMomMin = 5
ExtMomMax = 10       #*kF
ExtMomBin = 5       #0,1,2,...n-1
# qindex = [0, 1, 2, 3]
fig, ax = plt.subplots(3,2, sharex='col')

qindex = range(ExtMomBin)
# fig, ax = plt.subplots(3,3, sharex='col')

ax = ax.reshape(-1)

Para = param(D, Spin)
with open("./parameter", "r") as file:
    para = file.readline().split(" ")
    order = int(para[0])
    beta = float(para[1])
    rs   = float(para[2])
    # lam  = float(para[4])
fname = 'Data_Polar/polar_beta{0:3.1f}_rs{1:3.1f}_o{2}.dat'.format(beta, rs, order)
Polar = np.loadtxt(fname)
print('Loading '+fname)
fname = 'polar0_beta{0:.1f}_rs{1:.1f}.txt'.format(beta,rs)
Polar0 = np.loadtxt(fname)
print('Loading '+fname)

# print(Polar)

for index, qi in enumerate(qindex):
    # q = qi *Para.kF*1.5
    q = (index*(ExtMomMax-ExtMomMin)/ExtMomBin+ExtMomMin) *Para.kF
    print(index, q)

    polarization = {1:[], 2:[], 3:[], 4:[], 5:[], -1:[]}

    i=0
    for order in Polar[:,3]:
        if(abs(Polar[i,0]-q)<1e-4):
            polarization[int(order)].append(Polar[i])
        i = i+1
    i=0
    for q0 in Polar0[:,0]:
        if(abs(q0-q)<1e-4):
            polarization[-1].append(Polar0[i])
        i = i+1

    for o in range(1, Para.Order+1):
        polarization[o] = np.array(polarization[o])
        # print(o, polarization[o])

    # fig, ax = plt.subplots()
    ax1 = ax[index]
    if index in [0,1,2,3,4,5,6]:
        # print(polarization[-1][0][1])
        x = [0.45,3.05]
        y= [-polarization[-1][0][1],-polarization[-1][0][1]]
        for o in range(3, Para.Order+1):
            ax1.errorbar(polarization[o][:,4], polarization[o][:,1], yerr=polarization[o][:,2], fmt='o-', capthick=1, capsize=4,
                    color=ColorList[o], label="Order {0}".format(o))
        # ax1.errorbar(x, y, yerr=0, fmt='--', capthick=1, capsize=4, label=r"$\chi_0(q=q_0)$")
        ax1.set_xlabel(r"$\lambda/E_F$", size=size)
        ax1.legend(loc=1, frameon=False, fontsize=10)
    else:
        for o in range(3, Para.Order+1):
            ax1.errorbar(polarization[o][:,4], polarization[o][:,1], yerr=polarization[o][:,2], fmt='o-', capthick=1, capsize=4,
                    color=ColorList[o], label="Order {0}".format(o))
        ax1.legend(loc=1, frameon=False, fontsize=10)
    if index >= len(qindex)/2:
        ax1.set_xlabel(r"$\lambda/E_F$", size=size)
        # ax1.legend(loc=1, frameon=False, fontsize=size)
    if index in [0,2,4]:
        ax1.set_ylabel("$\Pi(q=q_0)/N_F$", size=size)

    ax1.set_title(r"$q_0={0:3.1f}k_F$".format(q/Para.kF))
    # plt.suptitle(r"$\beta={0}, r_s={1}, q_0={2:3.1f}k_F$".format(beta,rs,q/Para.kF))
    # plt.savefig('Plot/Polarvslam_charge_beta{0}_rs{1}_q{2:3.1f}.pdf'.format(beta,rs,q/Para.kF))

plt.savefig('Plot/Polarvslam_charge_beta{0}_rs{1}_0612.pdf'.format(beta,rs))
plt.show()


