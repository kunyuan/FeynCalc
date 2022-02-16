import numpy as np
import scipy.integrate as integrate

def __Ek0(k, Dim, Beta, Spin, Kf, mur):
    mu = mur*Kf*Kf
    x = Beta*(k*k-mu)
    f = Spin * k**2 /(1+np.exp(x))
    if Dim == 3:
        return f* k**2/2.0/np.pi**2
    elif Dim == 2:
        return f* k /2.0/np.pi
    else:
        return None


def Getmu(beta_):
    if beta_==0.0625:
        return -71.030575819501
    elif abs(beta_-0.08333333)<1e-6:
        return -48.067757106278
    elif beta_==0.125:
        return -27.136773591511
    elif abs(beta_-0.16666667)<1e-6:
        return -17.725431801201
    elif beta_ == 0.25:
        mu_r = -9.323691469984
    elif beta_ == 0.5:
        mu_r = -2.461438843431
    elif beta_ == 1:
        mu_r = -0.021460754987
    elif beta_ == 2:
        mu_r = 0.743112084259
    elif beta_ == 2.5:
        mu_r = 0.840347314915
    elif beta_ == 3:
        mu_r = 0.892499404946
    elif beta_ == 4:
        mu_r = 0.942615755201
    elif beta_ == 8:
        mu_r = 0.986801399943
    elif beta_ == 16:
        mu_r = 0.996768053583
    elif beta_ == 25:
        mu_r = 0.998680896718
    elif beta_ == 40:
        mu_r = 0.999485480206
    else:
        mu_r = 1
        print('beta>>1?  mu=Ef')    
    return mu_r

def Kinetic(Dim, Beta, Spin, Kf):
    assert Dim == 2 or Dim == 3, "Only Dim=2 and 3 are implemented."
    with open("./parameter", "r") as file:
        para = file.readline().split(" ")
        beta_ = float(para[1])
    mu_r = Getmu(beta_)

    return integrate.quad(__Ek0, 0.0, 100.0*Kf, args=(Dim, Beta, Spin, Kf, mu_r))