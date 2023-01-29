import numpy as np

###############################################################################
#### OUR ANALYTICAL RESULTS 1
###############################################################################
gradients     = np.array([0.10, 0.020, 0.005, 0.001])
gradients_unc =np.array([0.02, 0.005, 0.001, 0.001])

# Uncertainty propagation
coefsum = sum([(n+1)*gradients[n] for n in range(len(gradients))])
coefsum_unc = np.sqrt(
    sum([((n+1)*gradients_unc[n])**2 for n in range(len(gradients))]))
print(f"\nGradients weighted sum = {round(coefsum,5)} +- {round(coefsum_unc,5)}")

# Find probability distribution of n-lambdas
grad_sum = sum(gradients)
grad_sum_unc = np.sqrt(sum([g**2 for g in gradients_unc]))

lambd_probs = gradients/sum(gradients)
lambd_probs_uncs = np.sqrt(gradients_unc**2/grad_sum**2 +
                            gradients**2*grad_sum_unc**2/grad_sum**4)

P = sum([p*np.log(p) for p in lambd_probs])
P_unc = np.sqrt(sum((lambd_probs+1)**2 * lambd_probs_uncs**2))
C_hv = - sum(gradients)*sum([p*np.log(p) for p in lambd_probs])
C_hv_unc = np.sqrt(P**2 * grad_sum_unc**2 +
                    grad_sum**2 * P_unc**2)     
# Discreteness length in terms of Planckian length
l = 2*np.sqrt(C_hv)
l_unc = C_hv_unc/np.sqrt(C_hv)
print(f"\nDiscreteness scale as of ours = {round(l,5)} +- {round(l_unc,5)} l_p")


###############################################################################
#### OUR ANALYTICAL RESULTS 2
###############################################################################

# Find probability distribution of n-lambdas
grad_sum = sum(gradients)
grad_sum_unc = np.sqrt(sum([g**2 for g in gradients_unc]))

C_hv = sum(gradients*np.log(grad_sum/gradients))
dC_da_i = np.log(grad_sum/gradients)
C_hv_unc = np.sqrt( sum(gradients_unc**2 * dC_da_i**2))     
# Discreteness length in terms of Planckian length
l = 2*np.sqrt(C_hv)
l_unc = C_hv_unc/np.sqrt(C_hv)
print(f"\nA as of ours2 {grad_sum} +- {grad_sum_unc}")
print(f"C_hv as of ours2 {C_hv} +- {C_hv_unc}")
print(f"Discreteness scale as of ours2 = {round(l,5)} +- {round(l_unc,5)} l_p")


###############################################################################
#### UNCERTAINTY PACKAGE
###############################################################################
from uncertainties import ufloat
from uncertainties.umath import log, sqrt
from uncertainties import umath
from uncertainties import unumpy 


ai_s = unumpy.uarray(gradients, gradients_unc)
A = sum(ai_s)
C_hv = (ai_s * unumpy.log(A/ai_s)).sum()
l = 2*sqrt(C_hv)
print(f"\nA as of uncert {A}")
print(f"C_hv as of uncert {C_hv}")
print(f"Discreteness scale as of uncert = {l} l_p")

