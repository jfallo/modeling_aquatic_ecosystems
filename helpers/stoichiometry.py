import numpy as np
import pandas as pd

molarMass = {
    'C': 12,
    'O': 16,
    'H': 1,
    'N': 14,
    'P': 31,
    'S': 32
}

def mass_to_molar(e, mass):
    return mass / molarMass[e]

def molar_to_mass(e, mol):
    return mol * molarMass[e]

def get_stoich_coeffs_df(
        process, 
        compositions, 
        subs_and_orgs, 
        normalized_substance, 
        norm, 
        constraints, 
        additional_constraints= None):
    alpha = np.array([
        [compositions[s].get(c, 0.0) for s in subs_and_orgs]
        for c in constraints
    ])

    normalization_constraint = np.zeros(len(subs_and_orgs))
    normalization_constraint[subs_and_orgs.index(normalized_substance)] = norm

    rows = [alpha, normalization_constraint]
    if additional_constraints is not None:
        rows.extend(additional_constraints)

    A = np.vstack(rows)
    b = np.zeros(A.shape[0])
    b[len(constraints)] = 1.0

    coeffs = np.linalg.solve(A, b)

    return pd.DataFrame(
        [coeffs],
        index= [process],
        columns= subs_and_orgs
    )
