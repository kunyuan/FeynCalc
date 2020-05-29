import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat
import sys
import glob
import os
import re
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

rs = 1.0
Lambda = 1.0
Beta = 25

##############   3D    ##################################
# kF = (9.0*np.pi/4.0)**(1.0/3.0)/rs  # 3D
###### Bare Green's function    #########################
# Bubble=0.08871  # 3D, Beta=0.5, rs=1
# Bubble=0.0971916  #3D, Beta=10, rs=1
Bubble = 0.0971613  # 3D, T=0.04Ef, rs=1
# Bubble= 0.097226 # 3D, zero temperature, rs=1
###### Fock dressed Green's function ###################
# Bubble, Density = 0.088883, 0.2387  # 3D, Beta=0.1, rs=1, Lambda=1.0

##############   2D    ##################################
###### Bare Green's function    #########################
kF = np.sqrt(2.0)/rs  # 2D
# Bubble=0.11635  #2D, Beta=0.5, rs=1
Bubble = 0.15916  # 2D, Beta=10, rs=1

ScanOrder = [1, 2, 3]
Name = ["1", "2", "3", "4"]
# ScanOrder=[3]
Index = {}
Index[1] = []
Index[2] = []
Index[3] = []
Index[4] = []
# Index[5]=[1,]
DataAll = {}
Data = {}
DataOrderByOrder = {}
DataAtOrder = {}
Normalization = 1


folder = "./Beta{0}_rs{1}_lambda{2}/".format(Beta, rs, Lambda)

files = os.listdir(folder)
for order in ScanOrder:
    Num = 0
    data0 = None
    for f in files:
        if re.match("group"+Name[order-1]+"_pid[0-9]+.dat", f):
            print f
            Num += 1
            d = np.loadtxt(folder+f)
            if data0 is None:
                data0 = d
            else:
                data0[:, 1:] += d[:, 1:]

    print "Found {0} files.".format(Num)
    data0[:, 1:] /= Num

    # data0[:,1]*=(-1)**(order-1)

    # print data0

    DataAll[order] = np.array(data0)

Data[order] = []
for i in Index[order]:
    Num = 0
    data = None
    # for f in glob.glob("Diag"+str(order)+"_*_"+str(i)+".dat"):
    for f in files:
        if re.match("group"+str(order-1)+"_diag"+str(i)+"_pid[0-9]+.dat", f):
            print f
            Num += 1
            d = np.loadtxt(folder+f)
            # print f, d[0,1]
            if data is None:
                data = d
            else:
                data[:, 1:] += d[:, 1:]

        print "Found {0} files.".format(Num)
        data[:, 1:] /= Num
        # print data
        Data[order].append(np.array(data))

Normalization = DataAll[1][0, 1]/Bubble

for key in DataAll.keys():
    DataAll[key][:, 1] /= Normalization

# for i in range(len(Data[key])):
#     Data[key][i][:, 1] /= Normalization

for i in ScanOrder:
    DataOrderByOrder[i] = np.copy(DataAll[i])

# DataOrderByOrder[4][:,1]*=1.0/Density

for i in range(len(ScanOrder)):
    order = ScanOrder[i]
    DataAtOrder[order] = np.copy(DataOrderByOrder[ScanOrder[0]])
    for j in range(1, i):
        DataAtOrder[order][:, 1] += DataOrderByOrder[ScanOrder[j]][:, 1]


def ErrorPlot(p, d, color, marker, label=None, size=4, shift=False):
    data = np.array(d)
    data[:, 0] /= kF
    if shift:
        data[:, 1] -= data[0, 1]
    p.plot(data[:, 0], data[:, 1], marker=marker, c=color, label=label,
           lw=1, markeredgecolor="None", linestyle="--", markersize=size)
    # p.errorbar(data[:,0],data[:,1], yerr=data[:,2], c=color, ecolor=color, capsize=0, linestyle="None")
    # p.fill_between(data[:,0], data[:,1]-data[:,2], data[:,1]+data[:,2], alpha=0.5, facecolor=color, edgecolor=color)


w = 1-0.429

fig, ax = plt.subplots()
# ax=fig.add_axes()
# ax = fig.add_subplot(122)

# plt.subplot(1,2,2)
ColorList = ['k', 'r', 'b', 'g', 'm', 'c']

for i in range(0, len(ScanOrder)):
    o = ScanOrder[i]
    ErrorPlot(ax, DataOrderByOrder[o],
              ColorList[i], 's', "Order {0}".format(o))
    # ErrorPlot(ax, DataAtOrder[o], ColorList[i], 's', "Order {0}".format(o))

# ErrorPlot(ax, Data[1][0], 'k', 's', "Diag 1")
# ErrorPlot(ax, tmp, 'm', 's', "Diag 3+c 1")
# ErrorPlot(ax, DataAll[3], 'k', 'o', "Order 3")

# ErrorPlot(ax, Data[2][1], 'g', 'o', "Order 3 counterbubble 1")
# ErrorPlot(ax, Data[2][2], 'g', '*', "Order 3 counterbubble 2")
# ErrorPlot(ax, Data[2][3], 'g', '>', "Order 3 counterbubble 3")

# ErrorPlot(ax, Data[1][1], 'olive', 'o', "Order 3 shift 1")
# ErrorPlot(ax, Data[1][2], 'olive', '*', "Order 3 shift 2")

# ErrorPlot(ax, Data[3][0], 'k', 's', "Diag 1", shift=True)
# ErrorPlot(ax, Data[3][1], 'g', 's', "Diag 2", shift=True)
# ErrorPlot(ax, Data[3][2], 'r', '*', "Diag 3", shift=True)
# ErrorPlot(ax, Data[3][3], 'b', 's', "Diag 4", shift=True)
# ErrorPlot(ax, Data[3][4], 'olive', '*', "Diag 5", shift=True)
# ErrorPlot(ax, Data[3][5], 'm', 's', "Diag 6", shift=True)
# ErrorPlot(ax, Data[3][6], 'c', '*', "Diag 7", shift=True)

# ErrorPlot(ax, Data[5], 'g', 's', "Diag 6")


x = np.arange(0, 0.2, 0.001)
y = 0.5*x**w
# ax.plot(x,y,'k-', lw=2)

ax.set_xlim([0.0, DataAll[1][-1, 0]/kF])
# ax.set_xticks([0.0,0.04,0.08,0.12])
# ax.set_yticks([0.35,0.4,0.45,0.5])
# ax.set_ylim([-0.02, 0.125])
# ax.set_ylim([0.07, 0.125])
ax.set_xlabel("$q/k_F$", size=size)
# ax.xaxis.set_label_coords(0.97, -0.01)
# # ax.yaxis.set_label_coords(0.97, -0.01)
# ax.text(-0.012,0.52, "$-I$", fontsize=size)
ax.set_ylabel("$-P(\omega=0, q)$", size=size)

# ax.text(0.02,0.47, "$\\sim {\\frac{1}{2}-}\\frac{1}{2} {\\left( \\frac{r}{L} \\right)} ^{2-s}$", fontsize=28)

plt.legend(loc=1, frameon=False, fontsize=size)
# plt.title("2D density integral")
plt.tight_layout()

# plt.savefig("spin_rs1_lambda1.pdf")
plt.show()
