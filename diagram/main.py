from free_energy import *
from polar import *
import copy
import sys


def Generate(Order, VerOrder, SigmaOrder, IsSelfEnergy, IsSpinPolar, SPIN):
    LnZOrder = Order-1
    DiagFile = "./Diagram/HugenDiag{0}.txt".format(LnZOrder)
    LnZ = free_energy(LnZOrder)
    # Load pre-generated lnZ diagrams
    # build labeled Feynman diagram to unlabled Hugenholtz diagram mapping
    print "\nLoad Order {0} LnZ diagrams ...".format(LnZOrder)
    LnZ.LoadDiagrams(DiagFile)

    print red("\nThe optimimal LnZ diagrams:")
    OptLnZHugenDiagList = LnZ.OptimizeLoopBasis()

    Polar = polar(Order)

    UniqueUnLabelDiagList = []

    for d in OptLnZHugenDiagList:
        print "\n============================================================="
        print blue("Processing LnZ diagram: {0} with SymFactor: {1}".format(
            d.GetPermu(), d.SymFactor))

        print "Attach two external vertexes ..."
        OptPolarHugenDiagDict = Polar.AttachExtVer(d)

        # print "Check Tadpole..."
        # for p in OptPolarHugenDiagDict.keys():
        #     if diag.HasTadpole(p, Polar.GetReference()):
        #         del OptPolarHugenDiagDict[p]

        # print "Check Fock..."
        # for p in OptPolarHugenDiagDict.keys():
        #     if diag.HasFock(p, Polar.GetReference()):
        #         del OptPolarHugenDiagDict[p]

        print "Group polarization diagrams from the same LnZ diagram..."
        UnLabelDiagDeformList = Polar.Group(
            OptPolarHugenDiagDict, TimeRotation=True)
        # each element contains a deforamtion of hugenholz polarization
        # diagrams in the same LnZ group

        print red("Representative polarization Hugenholtz diagram:")
        for d in UnLabelDiagDeformList:
            diagram = copy.deepcopy(OptPolarHugenDiagDict[d[0]])
            diagram.SymFactor = diagram.SymFactor*len(d)
            UniqueUnLabelDiagList.append(diagram)
            print red("{0} with SymFactor {1}".format(
                diagram.GetPermu(), diagram.SymFactor))

    print yellow("Total Unique Polarization diagram: {0}".format(
        len(UniqueUnLabelDiagList)))

    print "Save diagrams ..."
    with open("./output/Diag{0}{1}_{2}_{3}.txt".format("Polar", Order, VerOrder, SigmaOrder), "w") as f:
        f.write(Polar.ToString(UniqueUnLabelDiagList,
                               VerOrder, SigmaOrder, IsSelfEnergy, IsSpinPolar, SPIN))


if __name__ == "__main__":
    # print "Input Diagram Order: "
    # Order = int(sys.argv[1])
    Order = 5
    IsSelfEnergy = False
    IsSpinPolar = True
    SPIN = 1
    for o in range(2, Order+1):
        for v in range(0, Order):
            for g in range(0, (Order-1)/2+1):
                if o+v+2*g > Order:
                    continue
                Generate(o, v, g, IsSelfEnergy, IsSpinPolar, SPIN)
