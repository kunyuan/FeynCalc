from IO import *
import scipy.integrate as integrate

D = 3
Spin = 2
Para = param(D, Spin)

# Data, Step, Groups, ReWeight,  Grids = LoadFile("data", "pid[0-9]+.dat")
# assert len(Groups) == Data[0].shape[0], "group num doesn't match the data num!"
# KGrid = Grids["KGrid"]
# ReWeight = np.array(ReWeight)
# print ReWeight
# # print Step

# Norm = [sum(data[0, :]) for data in Data]
# # Norm = [data[0, 0] for data in Data]

# Data = [np.sum(data[:, :], axis=1)/norm for (data, norm) in zip(Data, Norm)]
# Avg, Err = Estimate(Data, Step)
# print "Error:", Err

# weight = Err
# Z = sum(weight)
# ReWeight2 = np.array([r*w/Z for (r, w) in zip(ReWeight, weight)])
# ReWeight2[0] = ReWeight2[1]*2
# ReWeight = (ReWeight+ReWeight2)/2.0
# for (idx, g) in enumerate(Groups):
#     ReWeight[idx] = 4.0**g[0]

# ReWeight[0] = ReWeight[1]*2.0

# print Groups
# print ReWeight
# with open("reweight.data", "w") as f:
#     for r in ReWeight:
#         f.write("{0} ".format(r))

# Norm = [sum(data[0, :]) for data in Data]
# print Phys

# Accu = {}
# Each = {}
# for o in range(1, Para.Order+1):
#     Accu[o] = []
#     Each[o] = []

# for (d, data) in enumerate(Data):
#     polar = np.zeros([Para.Order, len(KGrid)])
#     for (idx, g) in enumerate(Groups):
#         if g == (0, ):
#             continue
#         order = g[0]+g[1]
#         if d == 1:
#             print "{0}_{1} => idx {2}, with order {3}".format(
#                 g[0], g[1], idx, order)
#         polar[order-1, :] += Data[d][idx, :]/Norm[d]

#     for o in range(1, Para.Order+1):
#         Each[o].append(polar[o-1, :])
#         Accu[o].append(np.sum(polar[0:o, :], axis=0))

# for o in Accu.keys():
#     Accu[o] = Estimate(Accu[o], Norm)
#     Each[o] = Estimate(Each[o], Norm)
