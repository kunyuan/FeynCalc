#!/usr/bin/python3
from IO import *
import reduce
import matplotlib as mat
import matplotlib.pyplot as plt

plt.switch_backend('TkAgg')
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

D = 3
Spin = 2

ExtMomMax = 5       #*kF
ExtMomBin = 10       #0,1,2,...n-1
qindex = range(1,ExtMomBin)

Para = param(D, Spin)
with open("./parameter", "r") as file:
    para = file.readline().split(" ")
    beta = float(para[1])
    rs   = float(para[2])
    # lam  = float(para[4])
Polar = np.loadtxt('polar.txt')

for index in qindex:
    q = index * ExtMomMax*Para.kF/ExtMomBin
    polarization = {1:[], 2:[], 3:[], 4:[], 5:[]}

    i=0
    for order in Polar[:,3]:
        if(abs(Polar[i,0]-q)<1e-5):
            polarization[int(order)].append(Polar[i])
        i = i+1

    for o in range(1, Para.Order+1):
        polarization[o] = np.array(polarization[o])
        print(o, polarization[o])

    fig, ax = plt.subplots()
    for o in range(2, Para.Order+1):
        plt.errorbar(polarization[o][:,4], polarization[o][:,1], yerr=polarization[o][:,2], fmt='o-', capthick=1, capsize=4,
                    color=ColorList[o], label="Order {0}".format(o))

    ax.set_ylabel("$-P(q=q_0)/N_F$", size=size)
    ax.set_xlabel(r"$\lambda/E_F$", size=size)
    plt.tight_layout()
    plt.legend(loc=3, frameon=False, fontsize=size)

    plt.suptitle(r"$\beta={0}, r_s={1}, q_0={2:3.1f}k_F$".format(beta,rs,q/Para.kF))
    plt.savefig('Plot/Polarvslam_charge_beta{0}_rs{1}_q{2:3.1f}.pdf'.format(beta,rs,q/Para.kF))
plt.show()
