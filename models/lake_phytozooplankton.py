import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from helpers.process_rate_factors import lim_Monod, rad_Monod_avg

params = {
    'ALG': {
        'k_growth': 0.5,    # 1/d
        'k_death': 0.1,     # 1/d
        'alpha_P': 0.003,   # gP/gDM
        'C_ini': 0.1,       # gDM/m3
        'beta': 0.046       # 1/degC
    },
    'ZOO': {
        'k_growth': 0.4,    # gDM/m3d
        'k_death': 0.05,    # 1/d
        'Y': 0.2,           # gDM/gDM
        'C_ini': 0.1,       # gDM/m3
        'beta': 0.08        # 1/degC
    },
    'HPO4': {
        'K': 0.002,         # gP/m3
        'C_in': 0.04,       # gP/m3
        'C_ini': 0.04       # gP/m3
    },
    'A': 5e+006,            # m2
    'h_epi': 5,             # m
    'Q_in': 5,              # m3/s
    'T': {
        'T0': 20,           # degC
        'min': 5,           # degC
        'max': 25           # degC
    },
    'I': {
        'K': 30,            # W/m2
        'min': 25,          # W/m2
        'max': 225          # W/m2
    },
    'lambda': {
        '1': 0.1,           # 1/m
        '2': 0.1            # m2/gDM
    },
    't_max': 230            # d
}

times = np.linspace(0, 4*365, 4*365 + 1)
y0 = [params['HPO4']['C_ini'], params['ALG']['C_ini'], params['ZOO']['C_ini']]


# --- Simple lake phyto- and zooplankton model with constant driving forces --- #

def rhs(y, t, par):
    C_HPO4, C_ALG, C_ZOO = y

    # system description
    V_ini = par['h_epi'] * par['A']
    flow = par['Q_in'] * 86400
    k_flush = flow / V_ini

    # Monod nutrient limitation
    f_N = lim_Monod(C_HPO4, par['HPO4']['K'])

    # process rates
    rho_growth_ALG = par['ALG']['k_growth'] * f_N * C_ALG
    rho_death_ALG = par['ALG']['k_death'] * C_ALG
    rho_growth_ZOO = par['ZOO']['k_growth'] * C_ALG * C_ZOO
    rho_death_ZOO = par['ZOO']['k_death'] * C_ZOO

    # transformation rates
    r_HPO4 = -par['ALG']['alpha_P'] * rho_growth_ALG
    r_ALG = rho_growth_ALG - rho_death_ALG - 1/par['ZOO']['Y'] * rho_growth_ZOO
    r_ZOO = rho_growth_ZOO - rho_death_ZOO

    # ODEs
    dC_HPO4 = k_flush * (par['HPO4']['C_in'] - C_HPO4) + r_HPO4
    dC_ALG = k_flush * -C_ALG + r_ALG
    dC_ZOO = k_flush * -C_ZOO + r_ZOO

    return [dC_HPO4, dC_ALG, dC_ZOO]


res_ode = odeint(rhs, y0, times, args= (params,))
C_HPO4 = res_ode[:, 0]
C_ALG = res_ode[:, 1]
C_ZOO = res_ode[:, 2]

# plot
plt.figure(figsize= (14,6))

plt.plot(times, C_HPO4, color= 'red', linestyle= '--', label= r'$C_{\mathrm{HPO4}}$')
plt.plot(times, C_ALG, color= 'black', linestyle= '-', label= r'$C_{\mathrm{ALG}}$')
plt.plot(times, C_ZOO, color= 'blue', linestyle= '--', label= r'$C_{\mathrm{ZOO}}$')
plt.ylim(0.0, 0.5)
plt.xlabel('Time [days]')
plt.ylabel(r'Concentration [gP/m³]')

plt.legend(loc= 'upper right')
plt.grid(True)
plt.tight_layout()
plt.show()


# --- Extension to periodic driving forces --- #

def rhs_ext(y, t, par):
    C_HPO4, C_ALG, C_ZOO = y
    
    # system description
    V_ini = par['h_epi'] * par['A']
    flow = par['Q_in'] * 86400
    k_flush = flow / V_ini

    # time dependent environmental conditions (light and temperature)
    omega = 2*np.pi / 365.25
    I0 = (
        0.5 * (par['I']['min'] + par['I']['max']) + 
        0.5 * (par['I']['max'] - par['I']['min']) * np.cos(omega * (t - par['t_max']))
    )
    T = (
        0.5 * (par['T']['min'] + par['T']['max']) + 
        0.5 * (par['T']['max'] - par['T']['min']) * np.cos(omega * (t - par['t_max']))
    )
    
    # Monod nutrient limitation
    lambda_ALG = par['lambda']['1'] + par['lambda']['2'] * C_ALG
    f_N_ALG = lim_Monod(C_HPO4, par['HPO4']['K'])
    f_T_ALG = np.exp(par['ALG']['beta'] * (T - par['T']['T0']))
    f_I_ALG = rad_Monod_avg(I0, par['I']['K'], lambda_ALG, par['h_epi'])

    f_T_ZOO = np.exp(par['ZOO']['beta'] * (T - par['T']['T0']))
    
    # process rates
    rho_growth_ALG = par['ALG']['k_growth'] * f_N_ALG * f_T_ALG * f_I_ALG * C_ALG
    rho_death_ALG = par['ALG']['k_death'] * C_ALG
    rho_growth_ZOO = par['ZOO']['k_growth'] * f_T_ZOO * C_ALG * C_ZOO
    rho_death_ZOO = par['ZOO']['k_death'] * C_ZOO
    
    # transformation rates
    r_HPO4 = -par['ALG']['alpha_P'] * rho_growth_ALG
    r_ALG = rho_growth_ALG - rho_death_ALG - 1/par['ZOO']['Y'] * rho_growth_ZOO
    r_ZOO = rho_growth_ZOO - rho_death_ZOO

    # ODEs
    dC_HPO4 = k_flush * (par['HPO4']['C_in'] - C_HPO4) + r_HPO4
    dC_ALG = k_flush * -C_ALG + r_ALG
    dC_ZOO = k_flush * -C_ZOO + r_ZOO

    return [dC_HPO4, dC_ALG, dC_ZOO]


params['ALG']['k_growth'] = 0.8

# 2 year simulation
res_ode_ext = odeint(rhs_ext, y0, times, args= (params,))
C_HPO4_ext = res_ode_ext[:, 0]
C_ALG_ext = res_ode_ext[:, 1]
C_ZOO_ext = res_ode_ext[:, 2]

# plot
plt.figure(figsize= (14,6))

plt.plot(times[:2*365 + 1], C_HPO4_ext[:2*365 + 1], color= 'red', linestyle= '--', label= r'$C_{\mathrm{HPO4}}$')
plt.plot(times[:2*365 + 1], C_ALG_ext[:2*365 + 1], color= 'black', linestyle= '-', label= r'$C_{\mathrm{ALG}}$')
plt.plot(times[:2*365 + 1], C_ZOO_ext[:2*365 + 1], color= 'blue', linestyle= '--', label= r'$C_{\mathrm{ZOO}}$')
plt.ylim(0.0, 1.9)
plt.xlabel('Time [days]')
plt.ylabel(r'Concentration [gP/m³]')

plt.legend(loc= 'upper right')
plt.grid(True)
plt.tight_layout()
plt.show()

# 1 year simulation plot
plt.figure(figsize= (14,6))

plt.plot(times[:365 + 1], C_HPO4_ext[:365 + 1], color= 'red', linestyle= '--', label= r'$C_{\mathrm{HPO4}}$')
plt.plot(times[:365 + 1], C_ALG_ext[:365 + 1], color= 'black', linestyle= '-', label= r'$C_{\mathrm{ALG}}$')
plt.plot(times[:365 + 1], C_ZOO_ext[:365 + 1], color= 'blue', linestyle= '--', label= r'$C_{\mathrm{ZOO}}$')
plt.ylim(0.0, 1.3)
plt.xlabel('Time [days]')
plt.ylabel(r'Concentration [gP/m³]')

plt.legend(loc= 'upper right')
plt.grid(True)
plt.tight_layout()
plt.show()
