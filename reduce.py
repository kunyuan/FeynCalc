import numpy as np


def Estimate(Data, Weights, axis=0):
    """ Return Mean and Error  with given weights"""
    # Assume weights are similar when calculating error bars
    Weights = np.array(Weights)
    Num = len(Weights)
    assert Num > 0, "Data is empty!"
    assert Data.shape[0] == Num, "Data and Weights size must match!"
    Avg = np.average(Data, weights=Weights, axis=0)
    Var = np.average((Data-Avg)**2, weights=Weights, axis=0)
    Err = np.sqrt(Var/(Num-1)) if Num > 1 else 0.0
    return np.array((Avg, Err))


def Reduce(Dict, Map):
    """reduce Dict.keys() to mapped keys"""
    mappedDict = {}
    for g in Dict.keys():
        if Map.has_key(g):
            key = Map[g]
            if mappedDict.has_key(key):
                mappedDict[key] += Dict[g]
            else:
                mappedDict[key] = Dict[g]
    return mappedDict


def EstimateGroup(DataDict, Steps, Phys, group):
    Norm = np.sum(DataDict[(0, )][:, :], axis=-1)  # shape=pid
    if DataDict.has_key(group):
        data = DataDict[group][:, :]/Norm[:, np.newaxis]*Phys
        return Estimate(data, Steps)
    else:
        return None


def GetData(Data, Groups, Steps, Phys, Map):
    Norm = np.sum(Data[:, 0, :], axis=-1)  # shape=pid
    Data = Data/Norm[:, np.newaxis, np.newaxis] * Phys
    # reduce (order, verorder, sigmaorder) to (order+verorder, sigmaorder)
    Dict = {}
    for (idx, g) in enumerate(Groups):
        if g == (0, ):
            continue
        Dict[g] = Data[:, idx, :]

    Dict = Reduce(Dict, Map)
    print "Groups {0} reduced to {1}".format(Groups, list(set(Dict.keys())))

    EsData = {}
    for key in sorted(Dict.keys()):
        data = Dict[key]
        # data = np.average(Dict[key], axis=-1)  # average external K
        y, err = Estimate(data, Steps, axis=0)
        EsData[key] = np.array((y, err))

    return EsData, Dict
