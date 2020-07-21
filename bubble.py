import numpy as np
import scipy.integrate as integrate


def __uniformbubble(e, Dim, Beta, Spin, Ef):
    x = Beta*(e-Ef)
    if Dim == 3:
        f = Beta*Spin/8.0/np.pi**2*e**0.5
    else:
        # Dim==2
        f = Beta*Spin/8.0/np.pi

    if x > -100.0 and x < 100.0:
        return f/(1.0+np.cosh(x))
    else:
        return f*np.exp(-abs(x))


def __bubble(k, q, Dim, Beta, Spin, Ef):
    x = Beta*(k*k-Ef)
    if Dim == 3:
        f = Spin/8.0/np.pi**2/q*k*np.log(((q**2-2.0*k*q)/(q**2+2.0*k*q))**2)
        if x > -100.0 and x < 100.0:
            return f/(1.0+np.exp(x))
        elif x < -100.0:
            return f
        else:
            return 0.0


def Bubble(Dim, Beta, Spin, Ef, Mom):
    assert Dim == 2 or Dim == 3, "Only Dim=2 and 3 are implemented."
    if(abs(Mom) < 1.0e-10):
        return integrate.quad(__uniformbubble, 0.0,
                              100.0*Ef, args=(Dim, Beta, Spin, Ef))
        # print "Polarization at Q={0}: {1}+-{2}".format(
        #     Mom, Bubble[0], Bubble[1])
    else:
        return integrate.quad(__bubble, 0.0, 100.0*Ef,
                              args=(Mom, Dim, Beta, Spin, Ef))
