import numpy as np
import pandas as pd
from helpers.convert import mass_to_molar, molar_to_mass


# --- simple stoichiometric solution --- #

compositions = {
    'HPO4': {
        'P': 1.0,   # gP/gHPO4-P
        'N': 0.0    # gN/gHPO4-P
    },
    'NH4': {
        'P': 0.0,   # gP/gN4-N
        'N': 1.0    # gN/gN4-N
    },
    'ALG': {
        'P': 0.01,  # gP/gALG
        'N': 0.1    # gN/gALG
    }
}

substances = list(compositions.keys())
elemental_constituents = list({                            
    e for comp in compositions.values()
    for e in comp
})


# define composition matrix
alpha = np.array([
    [compositions[s].get(e, 0.0) for s in substances]
    for e in elemental_constituents
])

alpha_df = pd.DataFrame(
    alpha,
    index= elemental_constituents,
    columns= substances
).round(3)


# fix a mass factor for normalization and solve for process stoichiometric coefficients
normalizer = np.zeros(len(substances))
normalizer[-1] = 1.0

A = np.vstack([alpha, normalizer])
b = normalizer

nu_growth = np.linalg.solve(A, b)
nu_resp = -nu_growth

# define stoichiometric coefficients matrix
nu = np.vstack([
    nu_growth,
    nu_resp
])

nu_df = pd.DataFrame(
    nu,
    index= ['growth of algae', 'respiration of algae'],
    columns= substances
).round(3)




# --- complex stoichiometric solution --- #

params = {
    'Y_ZOO': 0.2,  # gZOO/gALG
    'f_e': 0.4     # gPOM/gALG
}

elemental_constituents = ['H', 'N', 'O', 'P', 'C']
constraints = elemental_constituents + ['charge']

organism_compositions = {
    'ALG': {
        'O': 0.50,  # gO/gALG
        'H': 0.07,  # gH/gALG
        'N': 0.06,  # gN/gALG
        'P': 0.005  # gP/gALG
    },
    'ZOO': {
        'O': 0.50,  # gO/gZOO
        'H': 0.07,  # gH/gZOO
        'N': 0.06,  # gN/gZOO
        'P': 0.01   # gP/gZOO
    },
    'POM': {
        'O': 0.40,  # gO/gPOM
        'H': 0.07,  # gH/gPOM
        'N': 0.04,  # gN/gPOM
        'P': 0.007  # gP/gPOM
    }
}

organisms = list(organism_compositions.keys())

# choose carbon factors to guarantee factors sum to unity
for org in organisms:
    organism_compositions[org]['C'] = 1 - sum([
        organism_compositions[org][e]
        for e in elemental_constituents if e != 'C'
    ])

substance_compositions = {
    'NH4': {
        'N': 1,
        'H': 4,
        'charge': 1
    },
    'NO3': {
        'N': 1,
        'O': 3,
        'charge': -1
    },
    'HPO4': {
        'H': 1,
        'P': 1,
        'O': 4,
        'charge': -2
    },
    'HCO3': {
        'H': 1,
        'C': 1,
        'O': 3,
        'charge': -1
    },
    'O2': {
        'O': 2
    },
    'H': {
        'H': 1,
        'charge': 1
    },
    'H2O': {
        'H': 2,
        'O': 1
    }
}

substances = list(substance_compositions.keys())

# convert organisms from mass to molar
for org in organisms:
    organism_compositions[org] = {
        e: mass_to_molar(e, val)
        for e, val in organism_compositions[org].items()
    }

# compositions of substances and organisms combined
compositions = {**organism_compositions, **substance_compositions}
substances += organisms


# define composition matrix
alpha = np.array([
    [compositions[s].get(c, 0.0) for s in substances]
    for c in constraints
])

alpha_df = pd.DataFrame(
    alpha,
    index= constraints,
    columns= substances
).round(3)


def get_coeffs_df(process, subs_and_orgs, normalized_substance, norm, additional_constraints= None):
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


# - get stoichiometric coefficients for each process -
nu_df_rows = []

# algae growth with ammonium as nitrogen source
ALG_growth_NH4_substances = [s for s in substances if s not in {'NO3', 'ZOO', 'POM'}]
nu_df_rows.append(get_coeffs_df(
    'algae growth (NH4)',
    ALG_growth_NH4_substances,
    'ALG', 1.0
))

# algae growth with nitrate as nitrogen source
ALG_growth_NO3_substances = [s for s in substances if s not in {'NH4', 'ZOO', 'POM'}]
nu_df_rows.append(get_coeffs_df(
    'algae growth (NO3)',
    ALG_growth_NO3_substances,
    'ALG', 1.0
))

# algae respiration with ammonium as nitrogen source
nu_df_rows.append(get_coeffs_df(
    'algae respiration (NH4)',
    ALG_growth_NH4_substances,
    'ALG', -1.0
))

# algae respiration with nitrate as nitrogen source
nu_df_rows.append(get_coeffs_df(
    'algae respiration (NO3)',
    ALG_growth_NO3_substances,
    'ALG', -1.0
))

# algae death
ALG_death_substances = [s for s in substances if s not in {'NO3', 'ZOO'}]

Y_ALG_death = min(1, 
                  compositions['ALG']['N'] / compositions['POM']['N'],
                  compositions['ALG']['P'] / compositions['POM']['P'])

ALG_death_constraint = np.zeros(len(ALG_death_substances))
ALG_death_constraint[ALG_death_substances.index('POM')] = 1.0
ALG_death_constraint[ALG_death_substances.index('ALG')] = Y_ALG_death

nu_df_rows.append(get_coeffs_df(
    'algae death',
    ALG_death_substances,
    'ALG', -1.0,
    additional_constraints= [ALG_death_constraint]
))

# zooplankton growth with ammonium as nitrogen source
ZOO_growth_substances = [s for s in substances if s != 'NO3']

ZOO_growth_constraint = np.zeros(len(ZOO_growth_substances))
ZOO_growth_constraint[ZOO_growth_substances.index('ZOO')] = 1.0
ZOO_growth_constraint[ZOO_growth_substances.index('ALG')] = params['Y_ZOO']

POM_constraint = np.zeros(len(ZOO_growth_substances))
POM_constraint[ZOO_growth_substances.index('POM')] = 1.0
POM_constraint[ZOO_growth_substances.index('ALG')] = params['f_e']

nu_df_rows.append(get_coeffs_df(
    'zooplankton growth',
    ZOO_growth_substances,
    'ZOO', 1.0,
    additional_constraints= [ZOO_growth_constraint, POM_constraint]
))

# zooplankton respiration
ZOO_respiration_substances = [s for s in substances if s not in {'NO3', 'ALG', 'POM'}]
nu_df_rows.append(get_coeffs_df(
    'zooplankton respiration',
    ZOO_respiration_substances,
    'ZOO', -1.0
))

# zooplankton death
ZOO_death_substances = [s for s in substances if s not in {'NO3', 'ALG'}]

Y_ZOO_death = min(1, 
                  compositions['ZOO']['N'] / compositions['POM']['N'],
                  compositions['ZOO']['P'] / compositions['POM']['P'])

ZOO_death_constraint = np.zeros(len(ZOO_death_substances))
ZOO_death_constraint[ZOO_death_substances.index('POM')] = 1.0
ZOO_death_constraint[ZOO_death_substances.index('ZOO')] = Y_ZOO_death

nu_df_rows.append(get_coeffs_df(
    'zooplankton death',
    ZOO_death_substances,
    'ZOO', -1.0,
    additional_constraints= [ZOO_death_constraint]
))

nu_df = pd.concat(nu_df_rows, axis= 0).fillna(0).round(3).reindex(columns= substances)
print(nu_df)




# --- introducing sulfur through sulfate --- #

elemental_constituents += ['S']
constraints = elemental_constituents + ['charge']

# convert organisms back to mass
for org in organisms:
    organism_compositions[org] = {
        e: molar_to_mass(e, val)
        for e, val in organism_compositions[org].items()
    }


organism_compositions['ALG']['S'] = 0.003
organism_compositions['ZOO']['S'] = 0.003
organism_compositions['POM']['S'] = 0.003

# update carbon factors to guarantee factors sum to unity
for org in organisms:
    organism_compositions[org]['C'] = 1 - sum([
        organism_compositions[org][e]
        for e in elemental_constituents if e != 'C'
    ])

substance_compositions['SO4'] = {
    'S': 1,
    'O': 4,
    'charge': -2
}

# update substances list
substances = list(substance_compositions.keys())

# convert organisms from mass to molar
for org in organisms:
    organism_compositions[org] = {
        e: mass_to_molar(e, val)
        for e, val in organism_compositions[org].items()
    }

# compositions of substances and organisms combined
compositions = {**organism_compositions, **substance_compositions}
substances += organisms


# define composition matrix
alpha = np.array([
    [compositions[s].get(c, 0.0) for s in substances]
    for c in constraints
])

alpha_df = pd.DataFrame(
    alpha,
    index= constraints,
    columns= substances
).round(5)


# - get stoichiometric coefficients for each process -
nu_df_rows = []

# algae growth with ammonium as nitrogen source
ALG_growth_NH4_substances = [s for s in substances if s not in {'NO3', 'ZOO', 'POM'}]
nu_df_rows.append(get_coeffs_df(
    'algae growth (NH4)',
    ALG_growth_NH4_substances,
    'ALG', 1.0
))

# algae growth with nitrate as nitrogen source
ALG_growth_NO3_substances = [s for s in substances if s not in {'NH4', 'ZOO', 'POM'}]
nu_df_rows.append(get_coeffs_df(
    'algae growth (NO3)',
    ALG_growth_NO3_substances,
    'ALG', 1.0
))

# algae respiration with ammonium as nitrogen source
nu_df_rows.append(get_coeffs_df(
    'algae respiration (NH4)',
    ALG_growth_NH4_substances,
    'ALG', -1.0
))

# algae respiration with nitrate as nitrogen source
nu_df_rows.append(get_coeffs_df(
    'algae respiration (NO3)',
    ALG_growth_NO3_substances,
    'ALG', -1.0
))

# algae death
ALG_death_substances = [s for s in substances if s not in {'NO3', 'ZOO'}]

Y_ALG_death = min(1, 
                  compositions['ALG']['N'] / compositions['POM']['N'],
                  compositions['ALG']['P'] / compositions['POM']['P'],
                  compositions['ALG']['S'] / compositions['POM']['S'])

ALG_death_constraint = np.zeros(len(ALG_death_substances))
ALG_death_constraint[ALG_death_substances.index('POM')] = 1.0
ALG_death_constraint[ALG_death_substances.index('ALG')] = Y_ALG_death

nu_df_rows.append(get_coeffs_df(
    'algae death',
    ALG_death_substances,
    'ALG', -1.0,
    additional_constraints= [ALG_death_constraint]
))

# zooplankton growth with ammonium as nitrogen source
ZOO_growth_substances = [s for s in substances if s != 'NO3']

ZOO_growth_constraint = np.zeros(len(ZOO_growth_substances))
ZOO_growth_constraint[ZOO_growth_substances.index('ZOO')] = 1.0
ZOO_growth_constraint[ZOO_growth_substances.index('ALG')] = params['Y_ZOO']

POM_constraint = np.zeros(len(ZOO_growth_substances))
POM_constraint[ZOO_growth_substances.index('POM')] = 1.0
POM_constraint[ZOO_growth_substances.index('ALG')] = params['f_e']

nu_df_rows.append(get_coeffs_df(
    'zooplankton growth',
    ZOO_growth_substances,
    'ZOO', 1.0,
    additional_constraints= [ZOO_growth_constraint, POM_constraint]
))

# zooplankton respiration
ZOO_respiration_substances = [s for s in substances if s not in {'NO3', 'ALG', 'POM'}]
nu_df_rows.append(get_coeffs_df(
    'zooplankton respiration',
    ZOO_respiration_substances,
    'ZOO', -1.0
))

# zooplankton death
ZOO_death_substances = [s for s in substances if s not in {'NO3', 'ALG'}]

Y_ZOO_death = min(1, 
                  compositions['ZOO']['N'] / compositions['POM']['N'],
                  compositions['ZOO']['P'] / compositions['POM']['P'],
                  compositions['ZOO']['S'] / compositions['POM']['S'])

ZOO_death_constraint = np.zeros(len(ZOO_death_substances))
ZOO_death_constraint[ZOO_death_substances.index('POM')] = 1.0
ZOO_death_constraint[ZOO_death_substances.index('ZOO')] = Y_ZOO_death

nu_df_rows.append(get_coeffs_df(
    'zooplankton death',
    ZOO_death_substances,
    'ZOO', -1.0,
    additional_constraints= [ZOO_death_constraint]
))

nu_df = pd.concat(nu_df_rows, axis= 0).fillna(0).round(5).reindex(columns= substances)
print(nu_df)
