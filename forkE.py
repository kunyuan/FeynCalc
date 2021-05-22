#!/usr/bin/env python3
from scipy import interpolate
from scipy import integrate
from scipy import optimize
import numpy as np
import sys
import time


def freeE(k):
    return k**2/2/mass-EF

def v_q(q, lam):
    return 8*pi/(eps*(q**2+lam**2))

def Sigma_x(k):
    pp = 0
    if lam>0:
        pp = lam/kF*(np.arctan((k+kF)/lam)-np.arctan((k-kF)/lam))
    qq = 1 - pp - (lam**2+kF**2-k**2)/(4*k*kF)*np.log((lam**2+(k-kF)**2)/(lam**2+(k+kF)**2))
    return -2*kF/(pi*eps)*qq

def n_bose(x,beta):
    if x*beta>100:
        return 0.0
    elif x*beta<-100:
        return -1.0
    else:
        return 1./(np.exp(x*beta)-1.)
        
def ferm(x,beta):
    if x*beta>100:
        return 0.0
    elif x*beta<-100:
        return 1.0
    else:
        return 1./(np.exp(x*beta)+1.)

def ferm_0T(x):
    if x>0:
        return 0.0
    elif x<0:
        return 1.0

def Sigma0T_integrand(q):
    f1 = ferm_0T(freeE(q))
    f2 = np.log(((k-q)**2+lam**2)/((k+q)**2+lam**2))
    return q*f1*f2/pi/k

def Sigma0T(k):
    pp = 0
    if(lam>0):
        pp = lam/kF*(np.arctan((k+kF)/lam)-np.arctan((k-kF)/lam))
    qq = 1 - pp - (lam**2+kF**2-k**2)/(4*k*kF)*np.log((lam**2+(k-kF)**2)/(lam**2+(k+kF)**2))
    return -2*kF/(pi*eps)*qq

def Sigma_integrand(q, iter):
    if iter==0:
        f1 = ferm(q**2/2/mass - mu, beta)
    else:
        if q>kmax:
            f1 = ferm(q**2/2/mass + Sigma0T(q) - mu, beta)
        else:
            f1 = ferm(q**2/2/mass + sigma(q) - mu, beta)
    if (abs(k)<1e-10):
        f2 = -4*q/(q**2+lam**2)
        return q*f1*f2/pi
    else:
        f2 = np.log(((k-q)**2+lam**2)/((k+q)**2+lam**2))
        return q*f1*f2/pi/k


def ferm_density(x):
    tmp = np.array((km**2 + forkSigma - x)*beta, dtype=np.float64)
    ferm_integrand = km**2*np.exp(-tmp)/(np.exp(-tmp)+1.)/(2*np.pi**2)
    # ferm_integrand = km**2/(np.exp(tmp)+1.)/(2*np.pi**2)

    return integrate.romb(ferm_integrand, dx=kmax/MAXBIN) - 3/(8*pi)*rs**(-3)

def shift0_value(x):
    if x==0.125:
        return -27.136773591511
    elif x==0.25:
        return -9.323691469984
    elif x==0.5:
        return -2.461438843431
    elif x==1:
        return -0.021460754987
    elif x==2:
        return 0.743112084259
    elif x==4:
        return 0.942615755201
    elif x==8:
        return 0.986801399943
    elif x==16:
        return 0.996768053583
    else:
        str1 = input("beta>>1, set mu=EF ?(Y or N)")
        if str1 == 'Y' or str1 == 'y':
            return 1
        else:
            exit(-1)


if __name__ == '__main__':
    
    pi=np.pi

    beta0, rs, lam = (input("beta, rs, lam: ").split())
    beta0 = float(beta0)
    rs = float(rs)
    lam = float(lam)
    # beta0 = 16.0
    # rs = 1.0
    # lam = 3.0
    mass = 0.5
    eps = 1.
    kF = (9*pi/4.)**(1./3.) /rs
    print('kF=', kF)
    print('Beta=',beta0)
    print('Rs=',rs)
    print('Lambda=',lam)
    EF = kF**2

    shift0 = shift0_value(beta0)
    beta = beta0/EF
    mu = shift0*EF
    print('idel mu=', mu)

    kmax = 12*kF
    MAXBIN = 2**14  # must be 2**N
    iterNum = 7
    km=np.linspace(0, kmax, num=MAXBIN+1)
    # km=np.linspace(kmax/MAXBIN, kmax+kmax/MAXBIN, num=MAXBIN+1)

    size = 16
    forkSigma=[]
    fork = []
    Sigma_0T=[]
    sigma_order={}
    mu_iter = []

    forkSigma=np.zeros(MAXBIN+1)

    for i in range(iterNum):
        j = 0
        for k in np.nditer(km):
            # sigma = forkSigma[j]
            integral = integrate.quad(Sigma_integrand, 0, np.inf, args=(i,))
            forkSigma[j] = integral[0]
            j = j + 1
        sigma = interpolate.PchipInterpolator(km, forkSigma)
        mu0 = mu
        mu = optimize.fsolve(ferm_density, mu0)
        print('mu=',mu)  

#    filename = 'sigma/sigma_beta{0}_rs{1}.txt'.format(beta0,rs)
#    np.savetxt(filename, forkSigma, fmt='%.12e', delimiter=' ')
    filename = 'sigma/sigma3D_beta{0}_rs{1}_lam{2}.txt'.format(beta0,rs,lam)
    with open(filename,'w') as f:
        f.write(str(beta0)+" "+str(rs)+" "+str(lam)+"\n")
        f.write(str(kmax)+" "+str(MAXBIN+1)+"\n")
        f.write(str(mu[0])+" "+str(shift0)+"\n")
        np.savetxt(f, forkSigma,fmt='%.12e',newline=' ')

    # for k in np.nditer(km):
    #     integral = integrate.quad(Sigma0T_integrand, 0, np.inf)
    #     sigma = integral[0]
    #     Sigma_0T.append(sigma)
    # Sigma_0T = np.array(Sigma_0T)
    # for k in np.nditer(km):
    #     Sigma_0T.append(Sigma0T(k))
    # Sigma_0T = np.array(Sigma_0T)

    # shift0 = EF
    # shift0 = 0.998680896718   # beta=25
    # shift0 = 0.996768053583   # beta=16
    # shift0 = 0.994225922800   # beta=12
    # shift0 = 0.991641236370   # beta=10
    # shift0 = 0.986801399943   # beta=8
    # shift0 = 0.975963627735   # beta=6
    # shift0 = 0.942615755201   # beta=4
    # shift0 = 0.743112084259   # beta=2
    # shift0 = -0.021460754987  # beta=1
    # shift0 = -2.461438843431  # beta=0.5
    # shift0 = -9.323691469984  # beta=1/4
    # shift0 =-27.136773591511  # beta=1/8

