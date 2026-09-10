library("stoichcalc")


# --- simple stoichiometric solution --- #

HPO4 <- c(P = 1.0,   # gP/gHPO4-P
          N = 0.0)   # gN/gHPO4-P
NH4  <- c(P = 0.0,   # gP/gN4-N
          N = 1.0)   # gN/gN4-N
ALG  <- c(P = 0.01,  # gP/gALG
          N = 0.1)   # gN/gALG

subst.comp <- list(HPO4 = HPO4,
                   NH4  = NH4,
                   ALG  = ALG)

# build composition matrix
alpha <- calc.comp.matrix(subst.comp)
print(alpha)

# calculate stoichiometric coefficients
nu.growth.ALG <- calc.stoich.coef(alpha      = alpha,
                                  name       = "growth.ALG",
                                  subst      = c("HPO4", "NH4", "ALG"),
                                  subst.norm = "ALG",
                                  nu.norm    = 1)
nu.resp.ALG <- calc.stoich.coef(alpha      = alpha,
                                name       = "resp.ALG",
                                subst      = c("HPO4", "NH4", "ALG"),
                                subst.norm = "ALG",
                                nu.norm    = -1)
nu <- rbind(nu.growth.ALG, nu.resp.ALG)
print(nu)




# --- complex stoichiometric solution --- #

param <- list(alpha.O.ALG = 0.50,   # gO/gALG
              alpha.H.ALG = 0.07,   # gH/gALG
              alpha.N.ALG = 0.06,   # gN/gALG
              alpha.P.ALG = 0.005,  # gP/gALG
              alpha.O.ZOO = 0.50,   # gO/gZOO
              alpha.H.ZOO = 0.07,   # gH/gZOO
              alpha.N.ZOO = 0.06,   # gN/gZOO
              alpha.P.ZOO = 0.01,   # gP/gZOO
              alpha.O.POM = 0.40,   # gO/gPOM
              alpha.H.POM = 0.07,   # gH/gPOM
              alpha.N.POM = 0.04,   # gN/gPOM
              alpha.P.POM = 0.007,  # gP/gPOM
              Y.ZOO       = 0.2,    # gZOO/gALG
              f.e         = 0.4)    # gPOM/gALG

# define carbon factor such that factors sum to 1
param$alpha.C.ALG = 1 - (
    param$alpha.O.ALG
    + param$alpha.H.ALG
    + param$alpha.N.ALG
    + param$alpha.P.ALG
)
param$alpha.C.ZOO = 1 - (
    param$alpha.O.ZOO
    + param$alpha.H.ZOO
    + param$alpha.N.ZOO
    + param$alpha.P.ZOO
)
param$alpha.C.POM = 1 - (
    param$alpha.O.POM
    + param$alpha.H.POM
    + param$alpha.N.POM
    + param$alpha.P.POM
)

NH4  <- c(N = 1,
          H = 4,
          charge = 3)
NO3  <- c(N = 1,
          O = 3,
          charge = -1)
HPO4 <- c(H = 1,
          P = 1,
          O = 4,
          charge = -2)
HCO3 <- c(H = 1,
          C = 1,
          O = 3,
          charge = -1)
O2   <- c(O = 2,
          charge = 0)
H    <- c(H = 1,
          charge = 1)
H2O  <- c(H = 2,
          O = 1,
          charge = 0)

ALG <- c(O = param$alpha.O.ALG / 16,
         H = param$alpha.H.ALG / 1,
         N = param$alpha.N.ALG / 14,
         P = param$alpha.P.ALG / 31,
         C = param$alpha.C.ALG / 12)
ZOO <- c(O = param$alpha.O.ZOO / 16,
         H = param$alpha.H.ZOO / 1,
         N = param$alpha.N.ZOO / 14,
         P = param$alpha.P.ZOO / 31,
         C = param$alpha.C.ZOO / 12)
POM <- c(O = param$alpha.O.POM / 16,
         H = param$alpha.H.POM / 1,
         N = param$alpha.N.POM / 14,
         P = param$alpha.P.POM / 31,
         C = param$alpha.C.POM / 12)

subst.comp <- list(NH4  = NH4,
                   NO3  = NO3,
                   HPO4 = HPO4,
                   HCO3 = HCO3,
                   O2   = O2,
                   H    = H,
                   H2O  = H2O,
                   ALG  = ALG,
                   ZOO  = ZOO,
                   POM  = POM)

# build composition matrix
alpha.complex <- calc.comp.matrix(subst.comp)
print(alpha.complex)

# define additional params
param$Y.ALG.death = min(1,
                        param$alpha.N.ALG / param$alpha.N.POM,
                        param$alpha.P.ALG / param$alpha.P.POM)
param$Y.ZOO.death = min(1,
                        param$alpha.N.ZOO / param$alpha.N.POM,
                        param$alpha.P.ZOO / param$alpha.P.POM)

# calculate stoichiometric coefficients
nu.growth.ALG.NH4.complex <- calc.stoich.coef(alpha      = alpha.complex,
                                              name       = "growth.ALG.NH4",
                                              subst      = c("NH4", "HPO4", "HCO3", "O2", "H", "H2O", "ALG"),
                                              subst.norm = "ALG",
                                              nu.norm    = 1.0)
nu.growth.ALG.NO3.complex <- calc.stoich.coef(alpha      = alpha.complex,
                                              name       = "growth.ALG.NO3",
                                              subst      = c("NO3", "HPO4", "HCO3", "O2", "H", "H2O", "ALG"),
                                              subst.norm = "ALG",
                                              nu.norm    = 1.0)
nu.resp.ALG.complex <- calc.stoich.coef(alpha      = alpha.complex,
                                        name       = "resp.ALG",
                                        subst      = c("NH4", "HPO4", "HCO3", "O2", "H", "H2O", "ALG"),
                                        subst.norm = "ALG",
                                        nu.norm    = -1.0)
nu.death.ALG.complex <- calc.stoich.coef(alpha       = alpha.complex,
                                         name        = "death.ALG",
                                         subst       = c("NH4", "HPO4", "HCO3", "O2", "H", "H2O", "ALG", "POM"),
                                         subst.norm  = "ALG",
                                         nu.norm     = -1.0,
                                         constraints = list(c("ALG" = param$Y.ALG.death,
                                                              "POM" = 1.0)))
nu.growth.ZOO.complex <- calc.stoich.coef(alpha       = alpha.complex,
                                         name        = "growth.ZOO",
                                         subst       = c("NH4", "HPO4", "HCO3", "O2", "H", "H2O", "ALG", "ZOO", "POM"),
                                         subst.norm  = "ZOO",
                                         nu.norm     = 1.0,
                                         constraints = list(c("ZOO" = 1.0,
                                                              "ALG" = param$Y.ZOO),
                                                            c("POM" = 1.0,
                                                              "ALG" = param$f.e)))
nu.resp.ZOO.complex <- calc.stoich.coef(alpha       = alpha.complex,
                                        name        = "resp.ZOO",
                                        subst       = c("NH4", "HPO4", "HCO3", "O2", "H", "H2O", "ZOO"),
                                        subst.norm  = "ZOO",
                                        nu.norm     = -1.0)
nu.death.ZOO.complex <- calc.stoich.coef(alpha       = alpha.complex,
                                         name        = "death.ZOO",
                                         subst       = c("NH4", "HPO4", "HCO3", "O2", "H", "H2O", "ZOO", "POM"),
                                         subst.norm  = "ZOO",
                                         nu.norm     = -1.0,
                                         constraints = list(c("ZOO" = param$Y.ZOO.death,
                                                              "POM" = 1.0)))
nu.complex <- rbind(nu.growth.ALG.NH4.complex, 
                    nu.growth.ALG.NO3.complex,
                    nu.resp.ALG.complex,
                    nu.death.ALG.complex,
                    nu.growth.ZOO.complex,
                    nu.resp.ZOO.complex,
                    nu.death.ZOO.complex)
print(round(nu.complex, 3))




# --- introducing sulfur through sulfate --- #

param$alpha.S.ALG = 0.003
param$alpha.S.ZOO = 0.003
param$alpha.S.POM = 0.003

# update carbon factor such that factors sum to 1
param$alpha.C.ALG = 1 - (
    param$alpha.O.ALG
    + param$alpha.H.ALG
    + param$alpha.N.ALG
    + param$alpha.P.ALG
    + param$alpha.S.ALG
)
param$alpha.C.ZOO = 1 - (
    param$alpha.O.ZOO
    + param$alpha.H.ZOO
    + param$alpha.N.ZOO
    + param$alpha.P.ZOO
    + param$alpha.S.ZOO
)
param$alpha.C.POM = 1 - (
    param$alpha.O.POM
    + param$alpha.H.POM
    + param$alpha.N.POM
    + param$alpha.P.POM
    + param$alpha.S.POM
)

SO4 <- c(S      = 1,
         O      = 4,
         charge = -2)

ALG <- c(O = param$alpha.O.ALG / 16,
         H = param$alpha.H.ALG / 1,
         N = param$alpha.N.ALG / 14,
         P = param$alpha.P.ALG / 31,
         C = param$alpha.C.ALG / 12,
         S = param$alpha.S.ALG / 32)
ZOO <- c(O = param$alpha.O.ZOO / 16,
         H = param$alpha.H.ZOO / 1,
         N = param$alpha.N.ZOO / 14,
         P = param$alpha.P.ZOO / 31,
         C = param$alpha.C.ZOO / 12,
         S = param$alpha.S.ZOO / 32)
POM <- c(O = param$alpha.O.POM / 16,
         H = param$alpha.H.POM / 1,
         N = param$alpha.N.POM / 14,
         P = param$alpha.P.POM / 31,
         C = param$alpha.C.POM / 12,
         S = param$alpha.S.POM / 32)

subst.comp <- list(NH4  = NH4,
                   NO3  = NO3,
                   HPO4 = HPO4,
                   SO4  = SO4,
                   HCO3 = HCO3,
                   O2   = O2,
                   H    = H,
                   H2O  = H2O,
                   ALG  = ALG,
                   ZOO  = ZOO,
                   POM  = POM)

# build composition matrix
alpha.complex.SO4 <- calc.comp.matrix(subst.comp)
print(alpha.complex.SO4)

# update additional params
param$Y.ALG.death = min(1,
                        param$alpha.N.ALG / param$alpha.N.POM,
                        param$alpha.P.ALG / param$alpha.P.POM,
                        param$alpha.S.ALG / param$alpha.S.POM)
param$Y.ZOO.death = min(1,
                        param$alpha.N.ZOO / param$alpha.N.POM,
                        param$alpha.P.ZOO / param$alpha.P.POM,
                        param$alpha.S.ZOO / param$alpha.S.POM)

# calculate stoichiometric coefficients
nu.growth.ALG.NH4.complex.SO4 <- calc.stoich.coef(alpha      = alpha.complex.SO4,
                                                  name       = "growth.ALG.NH4",
                                                  subst      = c("NH4", "HPO4", "SO4", "HCO3", "O2", "H", "H2O", "ALG"),
                                                  subst.norm = "ALG",
                                                  nu.norm    = 1.0)
nu.growth.ALG.NO3.complex.SO4 <- calc.stoich.coef(alpha      = alpha.complex.SO4,
                                                  name       = "growth.ALG.NO3",
                                                  subst      = c("NO3", "HPO4", "SO4", "HCO3", "O2", "H", "H2O", "ALG"),
                                                  subst.norm = "ALG",
                                                  nu.norm    = 1.0)
nu.resp.ALG.complex.SO4 <- calc.stoich.coef(alpha      = alpha.complex.SO4,
                                            name       = "resp.ALG",
                                            subst      = c("NH4", "HPO4", "SO4", "HCO3", "O2", "H", "H2O", "ALG"),
                                            subst.norm = "ALG",
                                            nu.norm    = -1.0)
nu.death.ALG.complex.SO4 <- calc.stoich.coef(alpha       = alpha.complex.SO4,
                                             name        = "death.ALG",
                                             subst       = c("NH4", "HPO4", "SO4", "HCO3", "O2", "H", "H2O", "ALG", "POM"),
                                             subst.norm  = "ALG",
                                             nu.norm     = -1.0,
                                             constraints = list(c("ALG" = param$Y.ALG.death,
                                                                  "POM" = 1.0)))
nu.growth.ZOO.complex.SO4 <- calc.stoich.coef(alpha       = alpha.complex.SO4,
                                              name        = "growth.ZOO",
                                              subst       = c("NH4", "HPO4", "SO4", "HCO3", "O2", "H", "H2O", "ALG", "ZOO", "POM"),
                                              subst.norm  = "ZOO",
                                              nu.norm     = 1.0,
                                              constraints = list(c("ZOO" = 1.0,
                                                              "ALG" = param$Y.ZOO),
                                                            c("POM" = 1.0,
                                                              "ALG" = param$f.e)))
nu.resp.ZOO.complex.SO4 <- calc.stoich.coef(alpha       = alpha.complex.SO4,
                                            name        = "resp.ZOO",
                                            subst       = c("NH4", "HPO4", "SO4", "HCO3", "O2", "H", "H2O", "ZOO"),
                                            subst.norm  = "ZOO",
                                            nu.norm     = -1.0)
nu.death.ZOO.complex.SO4 <- calc.stoich.coef(alpha       = alpha.complex.SO4,
                                             name        = "death.ZOO",
                                             subst       = c("NH4", "HPO4", "SO4", "HCO3", "O2", "H", "H2O", "ZOO", "POM"),
                                             subst.norm  = "ZOO",
                                             nu.norm     = -1.0,
                                             constraints = list(c("ZOO" = param$Y.ZOO.death,
                                                                  "POM" = 1.0)))
nu.complex.SO4 <- rbind(nu.growth.ALG.NH4.complex.SO4, 
                        nu.growth.ALG.NO3.complex.SO4,
                        nu.resp.ALG.complex.SO4,
                        nu.death.ALG.complex.SO4,
                        nu.growth.ZOO.complex.SO4,
                        nu.resp.ZOO.complex.SO4,
                        nu.death.ZOO.complex.SO4)
print(signif(nu.complex.SO4, 3))