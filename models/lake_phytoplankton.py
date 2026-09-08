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
    'HPO4': {
        'K': 0.002,         # gP/m3
        'C_in': 0.04,       # gP/m3
        'C_ini': 0.004      # gP/m3
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
y0 = [params['HPO4']['C_ini'], params['ALG']['C_ini']]


# --- Simple lake phytoplankton model with constant driving forces --- #

def rhs(y, t, par):
    C_HPO4, C_ALG = y

    # system description
    V_ini = par['h_epi'] * par['A']
    flow = par['Q_in'] * 86400
    k_flush = flow / V_ini

    # Monod nutrient limitation
    f_N = lim_Monod(C_HPO4, par['HPO4']['K'])

    # process rates
    rho_growth = par['ALG']['k_growth'] * f_N * C_ALG
    rho_death = par['ALG']['k_death'] * C_ALG

    # transformation rates
    r_HPO4 = -par['ALG']['alpha_P'] * rho_growth
    r_ALG = rho_growth - rho_death

    # ODEs
    dC_HPO4 = k_flush * (par['HPO4']['C_in'] - C_HPO4) + r_HPO4
    dC_ALG = k_flush * -C_ALG + r_ALG

    return [dC_HPO4, dC_ALG]


res_ode = odeint(rhs, y0, times, args= (params,))
C_HPO4 = res_ode[:, 0]
C_ALG = res_ode[:, 1]

# plot
fig, ax = plt.subplots(2, 1, figsize= (8,7), sharex= True)

ax[0].plot(times, C_HPO4, color= 'blue')
ax[0].set_ylabel(r'$C_{\mathrm{HPO4}}$ [gP/m³]')
ax[0].grid(True)

ax[1].plot(times, C_ALG, color= 'green')
ax[1].set_ylabel(r'$C_{\mathrm{ALG}}$ [gDM/m³]')
ax[1].set_xlabel('Time [days]')
ax[1].grid(True)

plt.tight_layout()
plt.show()


# --- Extension to periodic driving forces --- #

def rhs_ext(y, t, par):
    C_HPO4, C_ALG = y

    V_epi = par['h_epi'] * par['A']
    flow = par['Q_in'] * 86400
    k_flush = flow / V_epi

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

    # limitation factors
    f_N = lim_Monod(C_HPO4, par['HPO4']['K'])                       # nutrient limitation
    f_T = np.exp(par['ALG']['beta'] * (T - par['T']['T0']))         # temperature limitation
    lambda_ = par['lambda']['1'] + par['lambda']['2'] * C_ALG
    f_I = rad_Monod_avg(I0, par['I']['K'], lambda_, par['h_epi'])   # avg light limitation

    # process rates
    rho_growth = par['ALG']['k_growth'] * f_N * f_T * f_I * C_ALG
    rho_death = par['ALG']['k_death'] * C_ALG

    # transformation rates
    r_HPO4 = -par['ALG']['alpha_P'] * rho_growth
    r_ALG = rho_growth - rho_death

    # ODEs
    dC_HPO4 = k_flush * (par['HPO4']['C_in'] - C_HPO4) + r_HPO4
    dC_ALG = k_flush * -C_ALG + r_ALG

    return [dC_HPO4, dC_ALG]


params['ALG']['k_growth'] = 0.8

# 4 year simulation
res_ode_ext = odeint(rhs_ext, y0, times, args= (params,))
C_HPO4_ext = res_ode_ext[:, 0]
C_ALG_ext = res_ode_ext[:, 1]

# plot
fig, ax = plt.subplots(2, 1, figsize= (8,7), sharex= True)

ax[0].plot(times, C_HPO4_ext, color= 'blue')
ax[0].set_ylabel(r'$C_{\mathrm{HPO4}}$ [gP/m³]')
ax[0].grid(True)

ax[1].plot(times, C_ALG_ext, color= 'green')
ax[1].set_ylabel(r'$C_{\mathrm{ALG}}$ [gDM/m³]')
ax[1].set_xlabel('Time [days]')
ax[1].grid(True)

plt.tight_layout()
plt.show()

# 1 year simulation plot
fig, ax = plt.subplots(2, 1, figsize= (8,7), sharex= True)

ax[0].plot(times[:365 + 1], C_HPO4_ext[:365 + 1], color= 'blue')
ax[0].set_ylabel(r'$C_{\mathrm{HPO4}}$ [gP/m³]')
ax[0].grid(True)

ax[1].plot(times[:365 + 1], C_ALG_ext[:365 + 1], color= 'green')
ax[1].set_ylabel(r'$C_{\mathrm{ALG}}$ [gDM/m³]')
ax[1].set_xlabel('Time [days]')
ax[1].grid(True)

plt.tight_layout()
plt.show()
