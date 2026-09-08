import numpy as np


# --- Temperature dependence factor --- #

def temp_exp(T, T0= 0, beta= 0.1):
    return np.exp(beta * (T - T0))


# --- Limitation by substance concentration --- #

def lim_Monod(C, K):
    return C / (K + C)

def lim_exp(C, K):
    return 1 - np.exp(-C/K)

def lim_Blackman(C, K):
    return np.minimum(C/K, 1)

def lim_Monodquad(C, K):
    return C**2 / (K**2 + C**2)


# --- Inhibition by substance concentration --- #

def inh_Monod(C, K):
    return K / (K + C)

def inh_exp(C, K):
    return np.exp(-C/K)

def inh_Blackman(C, K):
    return np.maximum(1 - C/K, 0)

def inh_Monodquad(C, K):
    return K**2 / (K**2 + C**2)


# --- Light dependence factor --- #

def rad_Monod(I, K):
    return I / (K + I)

def rad_Smith(I, K):
    return I / np.sqrt(K**2 + I**2)

def rad_Steele(I, I_opt):
    return I/I_opt * np.exp(1 - I/I_opt)


# --- Average light dependence factor (across depth) --- #

def rad_Monod_avg(I0, K, lambda_, h):
    return (
        1 / (lambda_ * h) 
        * np.log(
            (K + I0) /
            (K + I0 * np.exp(-lambda_ * h))
        )
    )

def rad_Smith_avg(I0, K, lambda_, h):
    return (
        1 / (lambda_ * h) 
        * np.log(
            (I0/K + np.sqrt(1 + (I0/K)**2)) /
            (I0 * np.exp(-lambda_ * h) / K - np.sqrt(1 + (I0 * np.exp(-lambda_ * h) / K)**2))
        )
    )

def rad_Steele_avg(I0, I_opt, lambda_, h):
    return (
        np.e / (lambda_ * h) 
        * (
            np.exp(- I0 * np.exp(-lambda_ * h) / I_opt) 
            - np.exp(- I0 / I_opt)
        )
    )
