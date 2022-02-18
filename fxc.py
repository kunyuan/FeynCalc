import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import numpy as np
plt.style.use(['science'])
plt.rcParams.update({'font.size': 18})

Nq = 4
oi = 3

q = [0.0, 0.1/Nq, 0.2/Nq, 0.3/Nq, 0.4/Nq]

data = {}
for qi in [2, 3, 4]:
    data[qi] = np.loadtxt(f"Plot/chi_beta100.0_rs1.0_lam4.0_o{oi}_q{qi}.dat1")

plt.figure(figsize=(8.0, 6.0))
style = ['o', 's', '*', 'v', '^']
# colors=['4477AA', 'EE6677', '228833', 'CCBB44', '66CCEE', 'AA3377', 'BBBBBB']
colors = ['k', 'r', 'b', 'g',  'y']
for i in [2, 3, 4, ]:
    d = data[i]
    # ratio=1.0-d[:,1]/d[:, 2]
    # err=(d[:,3]/d[:,2])*ratio
    # err=(d[:,3]/d[:,2]**2)*d[0,1]
    # err = d[:, 3]/d[:, 2]**2*d[0, 1]
    # err = gaussian_filter1d(err, 6)
    # plt.plot(d[:, 1], d[:, -1], style[i-1]+":", lw=1, markersize=3,
    #          color=colors[i-1], label=f"$order {i}$")
    plt.plot(d[:, 1], d[:, -1], style[i-1]+":", lw=1, markersize=3,
             color=colors[i-1], label=f"$q/k_F={q[i]}$")
    # plt.fill_between(d[:, 0], ratio-err, ratio+err,
    #                  alpha=0.1, edgecolor=colors[i-1], facecolor=colors[i-1],
    #                  linewidth=0, linestyle='dashdot', antialiased=True)

# ylist = [1.152, 1.296, 1.438, 1.576, 1.683]
# yerrlist = [0.002, 0.006, 0.009, 0.009, 0.015]
# for (yi, yd) in enumerate(ylist):
#     x = d[:, 0]
#     y = 1.0-1./np.array([yd for _ in x])
#     err = abs(yerrlist[yi]/yd**2)
#     yerr = np.array([err for _ in x])
#     plt.plot(x, y, markersize=1, color=colors[yi], alpha=0.6)
#     plt.fill_between(x, y-yerr, y+yerr,
#                      alpha=0.2, edgecolor=colors[yi], facecolor=colors[yi],
#                      linewidth=0, linestyle='dashdot', antialiased=True)


# plt.ylim([-10, 0.2])
plt.xlim([0.0, 0.4])
# plt.xticks(fontsize=10)
# plt.yticks(fontsize=10)
plt.legend(ncol=3)
plt.xlabel("$\\omega_n/E_F$")
# plt.text(1.8, 0.83, "VDiagMC")
# plt.ylabel("$1-\chi_s(q)/\chi_0(q)$")
# plt.ylabel("$G_-(q/k_F)/(q/q_{TF})^2$")
plt.ylabel("$f_{xc}^+(i\\omega_n)/N_F*(\\omega_n/q)^2$")
# plt.ylabel("$f_{xc}^{-}(q/k_F)\cdot (q_{TF}/8\pi)^2$")
plt.savefig("fxc_charge_rescale.pdf")
# plt.show()
