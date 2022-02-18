from color import *
import numpy as np
import glob
import re
import sys
import os

# import seaborn as sns

# import matplotlib as mat
# mat.use("agg")
# import matplotlib.pyplot as plt


# mat.rcParams.update({'font.size': 16})
# mat.rcParams["font.family"] = "Times New Roman"
# size = 12

# sns.set_style("whitegrid")
# sns.set_palette("colorblind", n_colors=16)


def GetLine(file):
    while True:
        line = file.readline().strip()
        if len(line) > 0 and line[0] != "#":
            return line


def getListOfFiles(dirName):
    listOfFiles = list()
    for (dirpath, dirnames, filenames) in os.walk(dirName):
        listOfFiles += [os.path.join(dirpath, file) for file in filenames]
    return listOfFiles


class param:
    # Order, Beta, Rs, Mass2, Lambda, Charge2, TotalStep = [None, ]*7
    # kF, Nf, EF, Bubble = [0.0, ]*4
    def __init__(self, D, Spin):
        # self.DataFolder = "Data"
        self.InputFile = "parameter"
        self.Dim = D
        self.Spin = Spin

        suffix = ["_eqTime", "_freq", "_freqTau", "_Ek"]

        with open(self.InputFile, "r") as file:
            para = file.readline().split(" ")
            self.Order = int(para[0])
            self.Beta = float(para[1])
            self.Rs = float(para[2])
            self.Mass2 = float(para[3])
            self.Lambda = float(para[4])
            self.MinExtMom = float(para[5])
            self.MaxExtMom = float(para[6])
            self.TotalStep = int(para[7])
            self.DataFolder = "beta{0}_rs{1}_lam{2}_o{3}".format(
                self.Beta, self.Rs, self.Lambda, self.Order) + suffix[int(para[-1])]
        print(self.DataFolder)

        if self.Dim == 3:
            self.kF = (9.0*np.pi/4.0)**(1.0/3.0)/self.Rs
            self.Nf = self.kF/4.0/np.pi**2*self.Spin
            # self.kF = (9.0*np.pi/4.0)**(1.0/3.0)/self.Rs
            # self.Nf = self.kF/2.0/np.pi**2*self.Spin
        elif self.Dim == 2:
            self.kF = np.sqrt(2.0)/self.Rs  # 2D
            self.Nf = 1.0/4.0/np.pi*self.Spin
        else:
            print("Not Implemented for Dimension {0}".format(self.Dim))
            sys.exit(0)

        #self.EF = self.kF**2/2.0
        self.EF = self.kF**2
        self.Beta /= self.EF
        self.MaxExtMom *= self.kF

        print(yellow("Parameters:"))
        print("Rs={0}, kF={1}, EF={2}, Beta={3}, Mass2={4}, Lambda={5}\n".format(
            self.Rs, self.kF, self.EF, self.Beta, self.Mass2, self.Lambda))

# For the given path, get the List of all files in the directory tree


def LoadFile(Folder, FileName):
    Groups = []
    ReWeight = []
    Step = []
    Data = []
    Grid = {}

    for f in getListOfFiles(Folder):
        if re.search(FileName, f) and not('group' in f) and not ('Diag' in f):
            # print ("Loading ", f)
            try:
                with open(f, "r") as file:
                    line = file.readline().strip().split(":")[1]
                    Step.append(float(line))
                    line = file.readline().strip().split(":")[1]
                    if len(Groups) == 0:
                        for e in [e for e in line.split(",") if len(e) > 0]:
                            Groups.append(tuple([int(o)
                                                 for o in e.split("_")]))

                    line = file.readline().strip().split(":")[1]
                    ReWeight = [float(e)
                                for e in line.split(",") if len(e) > 0]

                    while True:
                        g = file.readline().split(":")
                        if g[0].find("Grid") != -1:
                            key = g[0].strip(" #")
                            Grid[key] = np.fromstring(g[1], sep=' ')
                        else:
                            break

                data = np.loadtxt(f)
                # if abs(data[1,0]-2.948e8) > 5e6:
                # continue
                # print (f, data[1])
                Data.append(data)

            except Exception as e:
                print("Failed to load {0}".format(f))
                print(str(e))

    Data = np.array(Data)
    DataDict = {}
    for (idx, g) in enumerate(Groups):
        DataDict[g] = np.array(Data[:, idx, :])

    return DataDict, np.array(Step), Groups, np.array(ReWeight), Grid


def LoadFile_Diag(Folder):
    # Diags = {'0_0': None, '1_0_0_0': None, '1_0_1_0': None, '1_0_2_0': None, '1_0_2_1': None}
    DataDict = {}
    Step = {}
    Grid = None
    fname1 = "Diag"
    fname2 = "pid[0-9]+.dat"

    for f in getListOfFiles(Folder):
        if re.search(fname1, f) and re.search(fname2, f):
            # print ("Loading ", f)
            try:
                with open(f, "r") as file:
                    line = file.readline().strip().split(":")
                    DiagName = line[-2].split(",")[0]
                    DiagName = tuple([int(o)for o in DiagName.split("_")])
                    # Step.append(float(line[-1]))
                    data = np.loadtxt(file)
                if DiagName == (0, 0):
                    DiagName = (0,)
                if DiagName in DataDict:
                    # Diags[DiagName] = np.row_stack(Diags[DiagName],data[:,1])
                    DataDict[DiagName].append(data[:, 1])
                    Step[DiagName].append(float(line[-1]))
                else:
                    DataDict[DiagName] = [data[:, 1]]
                    Step[DiagName] = [float(line[-1])]
                Grid = data[:, 0]
            except Exception as e:
                print("Failed to load {0}".format(f))
                print(str(e))
    for g in DataDict.keys():
        DataDict[g] = np.array(DataDict[g])

    return DataDict, Step[0, ], Grid


def ErrorPlot(p, x, d, color='k', marker='s', label=None, size=4, shift=False):
    p.plot(x, d, marker=marker, c=color, label=label,
           lw=1, markeredgecolor="None", linestyle="--", markersize=size)


ColorList = ['k', 'r', 'b', 'g', 'm', 'c', 'navy',
             'y', 'cyan', 'darkgreen', 'violet', 'lime', 'purple']
ColorList = ColorList*40


if __name__ == '__main__':
    Para = param(3, 2)

    # dirName = "./data"
    dirName = Para.DataFolder
    filename = "pid[0-9]+.dat"

    LoadFile(dirName, filename)
    # for elem in getListOfFiles(dirName):
    #     print(elem)
