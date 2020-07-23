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


def __bubble3D(k, q, Dim, Beta, Spin, Kf):
    x = Beta*(k*k-Kf*Kf)
    f = Spin/8.0/np.pi**2/q*k*np.log(((q**2-2.0*k*q)/(q**2+2.0*k*q))**2)
    if x > -100.0 and x < 100.0:
        return f/(1.0+np.exp(x))
    elif x <= -100.0:
        return f
    else:
        return 0.0


def __bubble2D(k, q, Dim, Beta, Spin, Kf):
    """ Dim==2 """
    f = Spin/np.pi*k/q/np.sqrt(q**2-4.0*k**2)

    x = Beta*(k*k-Kf*Kf)
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


def Bubble(Dim, Beta, Spin, Kf, Mom):
    assert Dim == 2 or Dim == 3, "Only Dim=2 and 3 are implemented."
    if(abs(Mom) < 1.0e-10):
        return integrate.quad(__uniformbubble, 0.0,
                              100.0*Kf*Kf, args=(Dim, Beta, Spin, Kf*Kf))
        # print "Polarization at Q={0}: {1}+-{2}".format(
        #     Mom, Bubble[0], Bubble[1])
    else:
        if Dim == 3:
            return integrate.quad(__bubble3D, 0.0, 100.0*Kf,
                                  args=(Mom, Dim, Beta, Spin, Kf))
        else:
            # Dim==2
            result = integrate.quad(__bubble2D, 0.0, Mom/2.0,
                                    args=(Mom, Dim, Beta, Spin, Kf))
            # result = (result[0] + Spin/4.0/np.pi /
            #           (1.0+np.exp(-Beta*(Mom*Mom-Kf*Kf))), result[1])
            return result
