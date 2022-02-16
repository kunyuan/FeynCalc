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

def __uniformbubble(k, Dim, Beta, Spin, Kf, mur):
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


def __bubble3D(k, q, Dim, Beta, Spin, Kf, mur):
    mu = mur*Kf*Kf
    x = Beta*(k*k-mu)
    # x = Beta*(k*k-Kf*Kf)
    f = -Spin/8.0/np.pi**2/q*k*np.log(((q**2-2.0*k*q)/(q**2+2.0*k*q))**2)
    # Mass = 1.0
    # x = Beta*(k*k-Kf*Kf)/(2*Mass)
    # f = Spin*Mass/4.0/np.pi**2/q*k*np.log(((q**2-2.0*k*q)/(q**2+2.0*k*q))**2)
    # return f/(1.0+np.exp(x))
    if x > -200.0 and x < 200.0:
        return f/(1.0+np.exp(x))
    elif x <= -200.0:
        return f
    else:
        return 0.0


def __bubble2D(k, q, Dim, Beta, Spin, Kf, mur):
    """ Dim==2 """
    f = Spin/np.pi*k/q/np.sqrt(q**2-4.0*k**2)
    mu = mur *Kf*Kf
    x = Beta*(k*k-mu)
    if x > -100.0 and x < 100.0:
        return f/(1.0+np.exp(x))
    elif x <= -100.0:
        return f
    else:
        return 0.0

# def __bubble2D(k, q, Dim, Beta, Spin, Kf):
#     """ Dim==2 """
#     if k < q/2.0:
#         f = Spin*Beta/2.0/np.pi/q*k*np.sqrt(q**2-4.0*k**2)
#     else:
#         f = 0.0

#     x = Beta*(k*k-Kf*Kf)
#     if x > -50.0 and x < 50.0:
#         return f/(1.0+np.exp(x))/(1.0+np.exp(x))
#     elif x <= -50.0:
#         return f
#     else:
#         return 0.0

def __uniformbubbleCT(k, Dim, Beta, Spin, Kf, mur):
    mu = mur*Kf* Kf
    x = Beta*(k*k-mu)
    if Dim == 3:
        f = -Spin/4.0/np.pi**2
    else:
        # Dim==2
        f = -Spin/4.0/np.pi

    if x > -100.0 and x < 100.0:
        return f*(-Beta/2.0) / (1.0+np.cosh(x))
    else:
        return f*(-Beta/2.0) *np.exp(-abs(x))


def __bubbleCT3D(k, q, Dim, Beta, Spin, Kf, mur):
    mu = mur*Kf*Kf
    x = Beta*(k*k-mu)
    f = Spin/8.0/np.pi**2/q*k*np.log(((q**2-2.0*k*q)/(q**2+2.0*k*q))**2)

    if x > -200.0 and x < 200.0:
        return f*(-Beta/2.0)/(1.0+np.cosh(x))
    else:
        return f*(-Beta/2.0)*np.exp(-abs(x))

def __bubbleCT2D(k, q, Dim, Beta, Spin, Kf, mur):
    """ Dim==2 """
    f = Spin/np.pi*k/q/np.sqrt(q**2-4.0*k**2)
    mu = mur *Kf*Kf
    x = Beta*(k*k-mu)
    if x > -200.0 and x < 200.0:
        return f*(-Beta/2.0)/(1.0+np.cosh(x))
    else:
        return f*(-Beta/2.0)*np.exp(-abs(x))

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


def Bubble(Dim, Beta, Spin, Kf, Mom):
    assert Dim == 2 or Dim == 3, "Only Dim=2 and 3 are implemented."
    # bub = integrate.quad(__uniformbubble, 0.0,
    #                           100.0*Kf*Kf, args=(Dim, Beta, Spin, Kf*Kf))
    # print ("Polarization at Q={0}: {1}+-{2}".format(0, bub[0], bub[1]))     

    # mu_r = {[0.5,-2.461438843431], [1,-0.021460754987], [2,0.743112084259], [4,0.942615755201], [16,0.996768053583]} #[beta,mu*]
    with open("./parameter", "r") as file:
        para = file.readline().split(" ")
        beta_ = float(para[1])
    mu_r = Getmu(beta_)
    if(abs(Mom) < 1.0e-10):
        return integrate.quad(__uniformbubble, 0.0,
                            100.0*Kf, args=(Dim, Beta, Spin, Kf, mu_r))
        # return integrate.quad(__uniformbubble, 0.0,
                            #   100.0*Kf*Kf, args=(Dim, Beta, Spin, Kf*Kf, mu_r))
        # print "Polarization at Q={0}: {1}+-{2}".format(
        #     Mom, Bubble[0], Bubble[1])
    else:
        if Dim == 3:
            return integrate.quad(__bubble3D, 0.0, 100.0*Kf,
                                  args=(Mom, Dim, Beta, Spin, Kf, mu_r))
        else:
            # Dim==2
            result = integrate.quad(__bubble2D, 0.0, Mom/2.0,
                                    args=(Mom, Dim, Beta, Spin, Kf, mu_r))
            # result = (result[0] + Spin/4.0/np.pi /
            #           (1.0+np.exp(-Beta*(Mom*Mom-Kf*Kf))), result[1])
            return result

def BubbleCT(Dim, Beta, Spin, Kf, Mom):
    assert Dim == 2 or Dim == 3, "Only Dim=2 and 3 are implemented."
    with open("./parameter", "r") as file:
        para = file.readline().split(" ")
        beta_ = float(para[1])
    mu_r = Getmu(beta_)

    if(abs(Mom) < 1.0e-10):
        return integrate.quad(__uniformbubbleCT, 0.0,
                            100.0*Kf, args=(Dim, Beta, Spin, Kf, mu_r))
    else:
        if Dim == 3:
            return integrate.quad(__bubbleCT3D, 0.0, 100.0*Kf,
                                    args=(Mom, Dim, Beta, Spin, Kf, mu_r))
        else:
            # Dim==2
            result = integrate.quad(__bubbleCT2D, 0.0, Mom/2.0,
                                    args=(Mom, Dim, Beta, Spin, Kf, mu_r))
            # result = (result[0] + Spin/4.0/np.pi /
            #           (1.0+np.exp(-Beta*(Mom*Mom-Kf*Kf))), result[1])
            return result
