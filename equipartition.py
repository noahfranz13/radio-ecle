import pandas as pd
import numpy as np
from astropy import units as u, constants as c
from astropy.cosmology import Planck18 as cosmo


def _mini_converter(d_L, nu_p):
    '''
    Converts to the appropriate units for the paper
    '''
    return d_L/1e28, nu_p/1e10

def Gamma(F_p_mJy, nu_p, z, d_L, t_d, eta=1, f_V=0.36, f_A=1):
    """
    Barniol Duran Eq. 4
    """
    d_L, nu_p = _mini_converter(d_L, nu_p)
    prefactor = 12
    midfactor = F_p_mJy**(1/3) * d_L**(2/3) * nu_p**(-17/24) * eta**(35/72) * (1+z)**(-1/3) * t_d**(-17/24)
    postfactor = f_A**(-7/24) * f_V**(-1/24)
        
    return prefactor * midfactor * postfactor
    
def beta(t, R, z):
    """
    Barniol Duran Eq. 22
    """
    return ((c.c.cgs.value * t / (R*(1+z)) ) + 1)**-1

def Req(p, F_p_mJy, nu_p, z, d_L, epsilon_e=0.1, epsilon_B=None, gamma=1, f_V=0.36, f_A=1, nu_m=None):
    """
    Barniol Duran Eq. 27
    """
    
    # convert into proper units
    d_L28, nu_p10 = _mini_converter(d_L, nu_p)
    
    # deal with epsilon e and B corrections to the radius (end of sec. 3)
    # if epislon_B is None, assume equipartition model
    if epsilon_B is None:
        epsilon = 1 
    else:
        epsilon = 11/6 * (epsilon_B/epsilon_e)
    
    if gamma <= 1:
        epsilon_correction = epsilon**(1/17)
    else:
        # use the relativistic correction
        epsilon_correction = epsilon**(1/12)
        gamma = gamma*epsilon**(1/24)
        
    # hot proton correction factor (sec. 4.2.2)
    xi = 1 + epsilon_e**-1
    hot_proton_correction = xi**(1/(13+2*p))
    
    # correct for electrons that radiate at nu_m (sec. 4.2.1)
    # if nu_m is not None: # this still assumes nu_m < nu_a (almost always true for TDEs)
    #     y_e = gamma_e(F_p_mJy, d_L, nu_p, z, R_17, gamma=gamma, f_A=f_A)
    #     gamma_m = y_e * (nu_m/nu_p)**(1/2)
    # else:
    chi_e = (p-2)/(p-1) * epsilon_e * c.m_p.cgs.value/c.m_e.cgs.value
    gamma_m = np.max(np.array([
        np.ones(len(chi_e))*2 if isinstance(chi_e, (list, np.ndarray, pd.Series)) else 2,
        chi_e*(gamma-1)
    ]), axis=0)
    gamma_factor = gamma_m**((2-p)/(13+2*p))
    
    # if gamma is less than 1 we have to do a newtonian correction to these models
    # (beginning of sec. 3)
    newtonian_correction = 1
    if gamma <= 1:
        newtonian_correction = 4**(1/(13+2*p))
    
    # the traditional model without any of the above corrections
    prefactor = 1e17 * (21.8 * 525**(p-1))**(1/(13+2*p))
    midfactor = F_p_mJy**((6+p)/(13+2*p)) * d_L28**((2*(p+6))/(13+2*p)) * nu_p10**(-1) * (1+z)**-((19+3*p)/(13+2*p)) 
    postfactor= f_A**-((5+p)/(13+2*p)) * f_V**-(1/(13+2*p)) * gamma**((p+8)/(13+2*p))
    
    return prefactor * midfactor * postfactor * gamma_factor * newtonian_correction * epsilon_correction * hot_proton_correction

def Eeq(p, F_p_mJy, nu_p, z, d_L, epsilon_e=0.1, epsilon_B=None, gamma=1, f_V=0.36, f_A=1, nu_m=None):
    """
    Barniol Duran Eq. 28 
    """
    
    # convert into the proper units for Barniol-Duran
    d_L28, nu_p10 = _mini_converter(d_L, nu_p)
    
    # deal with epsilon e and B corrections to the energy (end of sec. 3)
    # if epislon_B is None, assume equipartition model 
    if epsilon_B is None:
        epsilon = 1
    else:
        epsilon = 11/6 * (epsilon_B/epsilon_e)
    
    if gamma <= 1:
        epsilon_correction = (11/17)*epsilon**(-6/17) + (6/17)*epsilon**(11/17)
    else:
        # use the relativistic correction
        epsilon_correction = (11/17)*epsilon**(-5/12) + (6/17)*epsilon**(7/12)
    
    # hot proton correction factor (sec. 4.2.2)
    xi = 1 + epsilon_e**-1
    hot_proton_correction = xi**(11/(13+2*p))
    
    # correct for electrons that radiate at nu_m (sec 4.2.1)
    # if nu_m is not None: # this still assumes nu_m < nu_a (almost always true for TDEs)
    #     y_e = gamma_e(F_p_mJy, d_L, nu_p, z, R_17, gamma=gamma, f_A=f_A)
    #     gamma_m = y_e * (nu_m/nu_p)**(1/2)
    # else:
    chi_e = (p-2)/(p-1) * epsilon_e * c.m_p.cgs.value/c.m_e.cgs.value
    gamma_m = np.max(np.array([
        np.ones(len(chi_e))*2 if isinstance(chi_e, (list, np.ndarray, pd.Series)) else 2,
        chi_e*(gamma-1)
    ]), axis=0)
    gamma_factor = gamma_m**(-11*(p-2)/(13+2*p))
    
    # correct the relativistic equations for newtonian 
    # (beginning of sec. 3)
    newtonian_correction = 1
    if gamma <= 1:
        newtonian_correction = 4**(11/(13+2*p))
    
    prefactor = 1.3e48 * 21.8**-(2*(p+1)/(13+2*p)) * 525**((11*(p-1)/(13+2*p)))
    midfactor = F_p_mJy**((14+3*p)/(13+2*p)) * d_L28**(2*(3*p+14)/(13+2*p)) * nu_p10**(-1) * (1+z)**(-(27+5*p)/(13+2*p))
    postfactor= f_A**(-3*(p+1)/(13+2*p)) * f_V**(2*(p+1)/(13+2*p)) * gamma**-((5*p+16)/(13+2*p))

    return prefactor * midfactor * postfactor * gamma_factor * newtonian_correction * epsilon_correction * hot_proton_correction

def gamma_e(F_p_mJy, d_L, nu_p, z, R_17, gamma=1, f_A=1, eta=1):
    """
    Barniol-Duran eq. 14
    """
    d_L, nu_p = _mini_converter(d_L, nu_p)
    
    prefactor = 525
    midfactor = F_p_mJy * d_L**2 * nu_p**-2 * eta**(5/3) * (1+z)**-3
    postfactor = gamma/(f_A*R_17**2)

    return prefactor * midfactor * postfactor
    
def N_e(p, F_p_mJy, d_L, nu_p, z, R_17, gamma=1, f_A=1, eta=1, epsilon_e=0.1, nu_m=None):
    """
    Assumes a spherical outflow, Eq. 15
    """
    
    y_e = gamma_e(F_p_mJy, d_L, nu_p, z, R_17, gamma=gamma, f_A=f_A, eta=eta)
    
    d_L, nu_p = _mini_converter(d_L, nu_p)
    
    chi_e = (p-2)/(p-1) * epsilon_e * c.m_p.cgs.value/c.m_e.cgs.value
    gamma_m = np.max(np.array([
        np.ones(len(chi_e))*2 if isinstance(chi_e, (list, np.ndarray, pd.Series)) else 2,
        chi_e*(gamma-1)
    ]), axis=0)

    # conversion factor from Cendes+2021 to isotropic N_e
    gamma_factor = 4 * (y_e/gamma_m)**(p-1)
    
    prefactor = 1e54
    midfactor = (F_p_mJy**3 * d_L**6 * nu_p**-5 * eta**(10/3) * (1+z)**-8)
    postfactor = 1/(f_A**2 * R_17**4)
        
    return gamma_factor * prefactor * midfactor * postfactor
    
    
def n_e(p, F_p_mJy, d_L, nu_p, z, R_17, gamma=1, f_A=1, eta=1, emitting_region_width_fraction=0.1, epsilon_e=0.1):
    
    # number of electrons
    N = N_e(p, F_p_mJy, d_L, nu_p, z, R_17, gamma=gamma, f_A=f_A, eta=eta, epsilon_e=0.1)

    # volume of emitting region
    R = 1e17 * R_17
    V = 4/3*np.pi * (R**3 - ( (1-emitting_region_width_fraction) * R )**3)
    
    # convert from N_e to n_e using the same method as Eftekhari+18
    # electron number density inside the shock
    n_e = N / V 
    
    # convert to electron number density outside the shock using the jump conditions
    n_ext = n_e / (4*gamma**2)
    
    return n_ext
