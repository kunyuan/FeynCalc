import seaborn as sns
import os
import sys
import re
import glob
import numpy as np
from color import *

import matplotlib.pyplot as plt
import matplotlib as mat
mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

sns.set_style("whitegrid")
sns.set_palette("colorblind", n_colors=16)


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
        self.DataFolder = "Data"
        self.InputFile = "inlist"
        self.Dim = D
        self.Spin = Spin

        with open(self.InputFile, "r") as file:
            para = file.readline().split(" ")
            self.Order = int(para[0])
            self.Beta = float(para[1])
            self.Rs = float(para[2])
            self.Mass2 = float(para[3])
            self.Lambda = float(para[4])
            self.MaxExtMom = float(para[5])
            self.TotalStep = int(para[6])

        if self.Dim == 3:
            self.kF = (9.0*np.pi/4.0)**(1.0/3.0)/self.Rs
            self.Nf = self.kF/4.0/np.pi**2*self.Spin
        elif self.Dim == 2:
            self.kF = np.sqrt(2.0)/self.Rs  # 2D
            print "Not Implemented for Dimension {0}".format(self.Dim)
            sys.exit(0)
        else:
            print "Not Implemented for Dimension {0}".format(self.Dim)
            sys.exit(0)

        self.EF = self.kF**2
        self.Beta /= self.EF
        self.MaxExtMom *= self.kF

        print yellow("Parameters:")
        print "Rs={0}, kF={1}, EF={2}, Beta={3}, Mass2={4}, Lambda={5}\n".format(
            self.Rs, self.kF, self.EF, self.Beta, self.Mass2, self.Lambda)

# For the given path, get the List of all files in the directory tree


def average(Data, Weights):
    Z = sum(Weights)
    return sum([d*w/Z for (d, w) in zip(Data, Weights)])


def Estimate(Data, Weights):
    """ Return Mean and Error  with given weights"""
    # Assume weights are similar when calculating error bars
    assert len(Data) == len(Weights), "Data and Weights size must match!"
    assert len(Weights) > 0, "Data is empty!"
    # Z = np.sum(Weights)
    Avg = average(Data, Weights)
    # Avg = sum(Data)/sum(Weights)
    if len(Data) > 1:
        # Var = sum((d/norm - Avg) ** 2*norm/Z for (d, norm)
        #           in zip(Data, Weights))
        Var = average([(d-Avg)**2 for d in Data], Weights)
        # Var = np.cov(Data, aweights=Weights)
        return Avg, np.sqrt(Var/(len(Data)-1))
    else:
        return Avg, Var*0.0


def LoadFile(Folder, FileName):
    Groups = []
    Norm = []
    Data = []
    Grid = {}

    for f in getListOfFiles(Folder):
        if re.search(FileName, f):
            print "Loading ", f
            try:
                with open(f, "r") as file:
                    line = file.readline().strip().split(":")[1]
                    for e in [e for e in line.split(",") if len(e) > 0]:
                        Groups.append(tuple([int(o) for o in e.split("_")]))
                    while True:
                        g = file.readline().split(":")
                        if g[0].find("Grid") != -1:
                            key = g[0].strip(" #")
                            Grid[key] = np.fromstring(g[1], sep=' ')
                        else:
                            break

                data = np.loadtxt(f)
                Data.append(data)

            except Exception as e:
                print "Failed to load {0}".format(f)
                print str(e)

    return Data, Groups, Grid


def ErrorPlot(p, x, d, color='k', marker='s', label=None, size=4, shift=False):
    p.plot(x, d, marker=marker, c=color, label=label,
           lw=1, markeredgecolor="None", linestyle="--", markersize=size)


ColorList = ['k', 'r', 'b', 'g', 'm', 'c', 'navy',
             'y', 'cyan', 'darkgreen', 'violet', 'lime', 'purple']
ColorList = ColorList*40


if __name__ == '__main__':
    Para = param(3, 2)

    dirName = "./data"
    filename = "pid[0-9]+.dat"

    LoadFile(dirName, filename)
    # for elem in getListOfFiles(dirName):
    #     print(elem)
