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