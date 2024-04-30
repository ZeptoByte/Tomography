# File That Contains Parameterization of the Atmospheric Flux


# Returns the Muon Flux in units of (m^(-2)s^(-1)sr^(-1)GeV/c(^-1))
# Return a flux for muons given an energy, height, and muon momentum

import numpy as np


# Temporary

#
def muon_flux_Gaisser(E_mu, theta):
    ep_pi = 115
    ep_k = 850
    
    return 1400*E_mu**(-2.7)*(1/(1+(1.1*E_mu*np.cos(theta)/ep_pi))+.054/(1+(1.1*E_mu*np.cos(theta)/ep_k)))

def muon_flux_Gaisser_standard(E_mu, cos_theta):
    ep_pi = 115
    ep_k = 850
    
    return 1400*E_mu**(-2.7)*(1/(1+(1.1*E_mu*cos_theta/ep_pi))+.054/(1+(1.1*E_mu*cos_theta/ep_k)))

# https://www.mpi-hd.mpg.de/gerda/inaug/poster/gerda_inauguration_muons.pdf
def muon_flux_approximate(E_mu, cos_theta):
    
    return 121.754*E_mu**(-2.7)*(cos_theta)**2



def muon_weights(E_mu, cos_theta, area, time):

    const = 150 # Integrated Flux per m2 per s
    weight_E =  E_mu**(-2.7)
    weight_theta = cos_theta**2

    weights = weight_E*weight_theta
    return weights*(150/np.sum(weights))*area*time

def gamma_from_KE(KE, m):
    gamma = (KE + m)/m
    return gamma

def beta_from_gamma(gamma):
    beta = np.sqrt(1-1/gamma**2)
    return beta
if __name__ == '__main__':
    # Test Code Here
    pass