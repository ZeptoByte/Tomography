# File that Contains Information about the Propagation of Muons through Material
# Two Components:
# Energy Loss (Ionization Losses) and Absorption (Muon Decays)
import numpy as np


# Bethe-Bloch-True
def bethe_bloch_precise():
    return 0
## Bethe-Bloch Approximate
# rho -> Density, B -> velocity, Z -> Charge of Paritcle/Proton Charge
def bethe_bloch_approx(p, B, Z):
    return p*2*Z**2/B**2

# Muon Range -> Until Zero Energy
def muon_range(Z, p, E, m):
    return (E-((m*E)/(m+E)))/(2*Z**2*p*2)

# Muon Range -> Distance Where the Muon is Still Above the detector Threshold Energy
def muon_range_detector(Z, p, E, E_f, m):
    return (E-((m*E)/(m+E)))/(2*Z**2*p*2) - (E_f-((m*E_f)/(m+E_f)))/(2*Z**2*p*2)

# Now, for computing radiation effects

def one_over_radiation_length(p, Znucl, Ar):

    alpha = 1/137
    r0 = 2.81794*10**(-15)*100 # m
    N = 6.022*10**23

    return 4*alpha*r0**2*p*N/Ar*Znucl*(1+Znucl)*np.log(183/(Znucl)**(1/3))

# L is length traversed through the material
def scattering_angle(L, Z, Znucl, P, B, p, Ar):
    c = 299792458*100
    #angle = (Z/(P*c*B))*20*np.sqrt(L*one_over_radiation_length(p, Znucl, Ar))

    # Momentum in MeV
    rms = (Z/(P*B))*20*np.sqrt(L*one_over_radiation_length(p, Znucl, Ar)*100)
    
    return rms


def compute_deflection_angle(rms):

    return np.random.normal(loc=0, scale=rms)



if __name__ == '__main__':
    pass