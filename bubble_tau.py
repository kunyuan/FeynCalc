## Polarization Pi(\tau) for external momentum q=0
import numpy as np
import scipy.integrate as integrate

# def __uniformbubble(e, Dim, Beta, Spin, Ef, mur):
#     mu = mur*Ef
#     x = Beta*(e-mu)
#     if Dim == 3:
#         f = Beta*Spin/8.0/np.pi**2*e**0.5
#     else:
#         # Dim==2
#         f = Beta*Spin/8.0/np.pi

#     if x > -200.0 and x < 200.0:
#         return f/(1.0+np.cosh(x))
#     else:
#         return f*np.exp(-abs(x))

#  Bubble with q=0, tau=0
def __normbubble(k, Dim, Beta, Spin, Kf, mur):
    mu = mur*Kf*Kf
    x = Beta*(k*k-mu)
    if Dim == 3:
        f = Spin/4.0/np.pi**2/Beta
    else:
        # Dim==2
        f = Spin/4.0/np.pi/Beta
    if x > -200.0 and x < 200.0:
        return f/(1.0+np.exp(x)) 
    elif x <= -200.0:
        return f
    else:
        return 0.0

def __uniformbubble(k, Dim, Beta, Spin, Kf, mur, extQ):
    mu = mur*Kf*Kf
    x = Beta*(k*k-mu)
    x1 = Beta*((k+extQ)**2-mu)
    x2 = Beta*((k-extQ)**2-mu)
    if Dim == 3:
        f = Spin/2.0/np.pi**2
        f1 = 1 + np.log((1+np.exp(-x1))/(1+np.exp(-x2)))/ (4*k*extQ*Beta)
    else:
        # Dim==2
        f = Spin/2.0/np.pi
    if x > -200.0 and x < 200.0:
        return f*k*k *f1 / (1.0+np.exp(x)) 
    else:
        return 0.0

def __uniformbubbleStatic(k, Dim, Beta, Spin, Kf, mur):
    mu = mur*Kf*Kf
    x = Beta*(k*k-mu)
    if Dim == 3:
        f = Spin/4.0/np.pi**2
    else:
        # Dim==2
        f = Spin/4.0/np.pi

    if x > -200.0 and x < 200.0:
        return f/(1.0+np.exp(x))
    elif x <= -200.0:
        return f
    else:
        return 0.0

def __bubble3DStatic(k, Dim, Beta, Spin, Kf, mur, q):
    mu = mur*Kf*Kf
    x = Beta*(k*k-mu)
    f = -Spin/8.0/np.pi**2/q*k*np.log(((q**2-2.0*k*q)/(q**2+2.0*k*q))**2)
    if x > -200.0 and x < 200.0:
        return f/(1.0+np.exp(x))
    elif x <= -200.0:
        return f
    else:
        return 0.0

def __bubble3D(k, Dim, Beta, Spin, Kf, mur, z, q):
    mu = mur*Kf*Kf
    x = Beta*(k*k-mu)
    f = -Spin/8.0/np.pi**2/q*k*np.log((z**2+(q**2-2.0*k*q)**2)/ (z**2+(q**2+2.0*k*q)**2))
    if x > -200.0 and x < 200.0:
        return f/(1.0+np.exp(x))
    elif x <= -200.0:
        return f
    else:
        return 0.0


def Getmu(beta_):
    if beta_ == 0.0625:
        mu_r = -71.030575819501
    elif beta_ == 0.125:
        mu_r = -27.136773591511
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


def Bubble(Dim, Beta, Spin, Kf, Tau, ExtQ):
    assert Dim == 2 or Dim == 3, "Only Dim=2 and 3 are implemented."
    # bub = integrate.quad(__uniformbubble, 0.0,
    #                           100.0*Kf*Kf, args=(Dim, Beta, Spin, Kf*Kf))
    # print ("Polarization at Q={0}: {1}+-{2}".format(0, bub[0], bub[1]))     

    # mu_r = {[0.5,-2.461438843431], [1,-0.021460754987], [2,0.743112084259], [4,0.942615755201], [16,0.996768053583]} #[beta,mu*]
    with open("./parameter", "r") as file:
        para = file.readline().split(" ")
        beta_ = float(para[1])
    mu_r = Getmu(beta_)

    if (abs(Tau) < 1e-10 and abs(ExtQ) < 1e-10):
        return integrate.quad(__normbubble, 0.0,
                    100.0*Kf, args=(Dim, Beta, Spin, Kf, mu_r))
    else:
        return integrate.quad(__uniformbubble, 0.0,
                        100.0*Kf, args=(Dim, Beta, Spin, Kf, mu_r, ExtQ))

def Bubble_wn(Dim, Beta, Spin, Kf, OmegaN, ExtQ):
    assert Dim == 2 or Dim == 3, "Only Dim=2 and 3 are implemented."
    # bub = integrate.quad(__uniformbubble, 0.0,
    #                           100.0*Kf*Kf, args=(Dim, Beta, Spin, Kf*Kf))
    # print ("Polarization at Q={0}: {1}+-{2}".format(0, bub[0], bub[1]))     

    # mu_r = {[0.5,-2.461438843431], [1,-0.021460754987], [2,0.743112084259], [4,0.942615755201], [16,0.996768053583]} #[beta,mu*]
    with open("./parameter", "r") as file:
        para = file.readline().split(" ")
        beta_ = float(para[1])
    mu_r = Getmu(beta_)

    if (abs(OmegaN) < 1e-10 and abs(ExtQ) < 1e-10):
        return integrate.quad(__uniformbubbleStatic, 0.0,
                    100.0*Kf, args=(Dim, Beta, Spin, Kf, mu_r))
    elif (abs(OmegaN) < 1e-10):
        return integrate.quad(__bubble3DStatic, 0.0,
                        100.0*Kf, args=(Dim, Beta, Spin, Kf, mu_r, ExtQ))
        # return integrate.quad(__bubble3D, 0.0,
                        # 100.0*Kf, args=(Dim, Beta, Spin, Kf, mu_r, ExtQ, OmegaN), points=[-ExtQ/2,ExtQ/2])
    else:
        return integrate.quad(__bubble3D, 0.0,
                        200.0*Kf, args=(Dim, Beta, Spin, Kf, mu_r, OmegaN, ExtQ))
