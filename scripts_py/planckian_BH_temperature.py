import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpmath as mpm

params = {'text.usetex' : True,
          'font.size' : 20,
          'font.family' : 'Times New Roman',
          'axes.labelsize':15,
          'legend.fontsize': 20,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15,
          'figure.figsize': [8.5, 6.5], 
          'axes.prop_cycle':plt.cycler(color=
                            plt.rcParams['axes.prop_cycle'].by_key()['color']
                            +['magenta'])
          }
plt.rcParams.update(params)

def T_inv_zetaf(mu, ell = 1, varsigma1 = 1, x = 2, a=-0.1):
    if not hasattr(mu, "__len__"):
        return 8*np.pi*ell/mu + varsigma1 + mu * np.abs(a)*mpm.polylog(x-1, mu)
    else:
        res = []
        for mu_i in mu:
            res.append(8*np.pi*ell/mu_i + varsigma1 + mu_i * np.abs(a)*mpm.polylog(x-1, mu_i))
        return np.array(res)

ell = 1
varsigma_1 = 1
varsigma_2 = 0.1
a = -0.1
x = 1.5
mu_d = 1


M = np.logspace(np.log10(1/mu_d)+0.00001, 38, 100)
mu = 1/M
plt.figure(tight_layout = True)
plt.plot(1/mu, 1/T_inv_zetaf(mu, ell, varsigma_1, x, a), '.', ls = "")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$M \; [M_P]$")
plt.ylabel(r"$T \; [T_P]$")


M = np.logspace(np.log10(1/mu_d)+0.001, 1, 1000)
mu = 1/M
plt.figure(tight_layout = True)
plt.plot(1/mu, 1/T_inv_zetaf(mu, ell, varsigma_1, x, a), '.', ls = "")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$M \; [M_P]$")
plt.ylabel(r"$T \; [T_P]$")


M = np.logspace(np.log10(1/mu_d)+0.0000001, 0.00009, 1000)
mu = 1/M
plt.figure(tight_layout = True)
plt.plot(1/mu, 1/T_inv_zetaf(mu, ell, varsigma_1, x, a), '.', ls = "")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$M \; [M_P]$")
plt.ylabel(r"$T \; [T_P]$")
plt.show()

# plt.figure(tight_layout = True)
# plt.plot(mu, 1/T_inv_zetaf(mu, ell, varsigma_1, x, a), '.', ls = "")
# plt.yscale("log")
# plt.xscale("log")
# plt.xlabel(r"$\mu [\ell]$")
# plt.ylabel(r"Temperature")
# plt.show()




