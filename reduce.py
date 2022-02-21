import numpy as np
import math

def Estimate(Data, Weights, axis=0):
    """ Return Mean and Error  with given weights"""
    # Assume weights are similar when calculating error bars
    Weights = np.array(Weights)
    Num = len(Weights)
    assert Num > 0, "Data is empty!"
    assert Data.shape[0] == Num, "Data and Weights size must match!"
    Avg = np.average(Data, weights=Weights, axis=0)
    Var = np.average((Data-Avg)**2, weights=Weights, axis=0)
    Err = np.sqrt(Var/(Num-1)) if Num > 1 else np.zeros(len(Avg))
    return np.array((Avg, Err))


def Reduce(Dict, Map):
    """reduce Dict.keys() to mapped keys"""
    mappedDict = {}
    for g in Dict.keys():
        if g in Map:
            key = Map[g]
            if key in mappedDict:
                mappedDict[key] += Dict[g]
            else:
                mappedDict[key] = Dict[g]
    return mappedDict


def EstimateGroup(DataDict, Steps, Phys, group):
    Norm = np.sum(DataDict[(0, )][:, :], axis=-1)  # shape=pid
    if group in DataDict:
        data = DataDict[group][:, :]/Norm[:, np.newaxis]*Phys
        Newdata = []
        NewSteps = []
        for (idx, g) in enumerate(data[:,0]):
            # print(idx, data[idx, 0])
            if True in np.isnan(data[idx, :]):
                print("Here Nan!", data[idx,:])
                # np.delete(data, idx, axis=0)
            elif True in np.isinf(data[idx, :]):
                print("Here Inf!", data[idx, :])
            else:
                Newdata.append(data[idx,:])
                NewSteps.append(Steps[idx])
        Newdata = np.array(Newdata)
        NewSteps = np.array(NewSteps)
        # if group == (1, 0, 0, 2):
        #     print(group, Newdata)
        # print(group, Newdata)
        return Estimate(Newdata, NewSteps)
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
    print ("Groups {0} reduced to {1}".format(Groups, list(set(Dict.keys()))))

    EsData = {}
    for key in sorted(Dict.keys()):
        data = Dict[key]
        # data = np.average(Dict[key], axis=-1)  # average external K
        y, err = Estimate(data, Steps, axis=0)
        EsData[key] = np.array((y, err))

    return EsData, Dict
