library("ecosim")
library("stoichcalc")


# --- model with constant driving forces --- #

# set parameters
param <- list(alpha.O.ALG          = 0.50,    # gO/gALG
              alpha.H.ALG          = 0.07,    # gH/gALG
              alpha.N.ALG          = 0.06,    # gN/gALG
              alpha.P.ALG          = 0.005,   # gP/gALG
              alpha.O.ZOO          = 0.50,    # gO/gZOO
              alpha.H.ZOO          = 0.07,    # gH/gZOO
              alpha.N.ZOO          = 0.06,    # gN/gZOO
              alpha.P.ZOO          = 0.01,    # gP/gZOO
              alpha.O.POM          = 0.39,    # gO/gPOM
              alpha.H.POM          = 0.07,    # gH/gPOM
              alpha.N.POM          = 0.06,    # gN/gPOM
              alpha.P.POM          = 0.007,   # gP/gPOM
              Y.ZOO                = 0.2,     # gZOO/gALG
              f.e                  = 0.2,     # gPOM/gALG
              f.I                  = 0.2,     # gPOMI/(gPOMD + gPOMI)
              k.growth.ALG         = 0.8,     # 1/d
              k.growth.ZOO         = 0.4,     # m3/gDM/d
              k.resp.ALG           = 0.1,     # 1/d
              k.resp.ZOO           = 0.1,     # 1/d
              k.death.ALG          = 0.1,     # 1/d
              k.death.ZOO          = 0.05,    # 1/d
              k.nitri              = 0.1,     # gN/m3/d
              k.miner.ox.POM       = 0.02,    # 1/d
              k.miner.ox.POM.sed   = 5.0,     # gDM/m2/d
              k.miner.anox.POM.sed = 5.0,     # gDM/m2/d
              K.HPO4               = 0.002,   # gP/m3
              K.N                  = 0.04,    # gN/m3
              p.NH4                = 5,       # -
              K.O2.ZOO             = 0.2,     # gO/m3
              K.NH4.nitri          = 0.5,     # gN/m3
              K.O2.resp            = 0.5,     # gO/m3
              K.O2.nitri           = 0.4,     # gO/m3
              K.O2.miner           = 0.5,     # gO/m3
              K.POM.miner.sed      = 10,      # gDM/m2
              K.NO3.miner          = 0.1,     # gN/m3
              A                    = 5e+006,  # m2
              h.epi                = 5,       # m
              h.hypo               = 10,      # m
              Q.in                 = 5,       # m3/s
              C.NO3.in             = 0.5,     # gN/m3
              C.HPO4.in            = 0.04,    # gP/m3
              C.O2.in              = 10,      # gO/m3
              C.ALG.ini            = 0.1,     # gDM/m3
              C.ZOO.ini            = 0.1,     # gDM/m3
              C.POMD.ini           = 0.0,     # gDM/m3
              C.POMI.ini           = 0.0,     # gDM/m3
              C.NH4.ini            = 0.1,     # gN/m3
              C.NO3.ini            = 0.5,     # gN/m3
              C.HPO4.ini           = 0.04,    # gP/m3
              C.O2.ini             = 10,      # gO/m3
              D.POMD.ini           = 0.0,     # gDM/m2
              D.POMI.ini           = 0.0,     # gDM/m2
              beta.ALG             = 0.046,   # 1/degC
              beta.ZOO             = 0.08,    # 1/deg
              beta.BAC             = 0.046,   # 1/deg
              T0                   = 20,      # degC
              T.min                = 5,       # degC
              T.max                = 25,      # degC
              K.I                  = 30,      # W/m2
              I0.min               = 25,      # W/m2
              I0.max               = 225,     # W/m2
              lambda.1             = 0.10,    # 1/m
              lambda.2             = 0.10,    # m2/gDM
              v.ex.O2              = 1,       # m/d
              v.sed.POM            = 1,       # m/d
              Kz.summer            = 0.02,    # m2/d
              Kz.winter            = 20,      # m2/d
              h.meta               = 5,       # m
              t.max                = 230,     # d
              p                    = 101325)  # Pa

param$alpha.C.ALG <- 1 - (param$alpha.O.ALG
                           + param$alpha.H.ALG
                           + param$alpha.N.ALG
                           + param$alpha.P.ALG)
param$alpha.C.ZOO <- 1 - (param$alpha.O.ZOO
                           + param$alpha.H.ZOO
                           + param$alpha.N.ZOO
                           + param$alpha.P.ZOO)
param$alpha.C.POM <- 1 - (param$alpha.O.POM
                           + param$alpha.H.POM
                           + param$alpha.N.POM
                           + param$alpha.P.POM)

param$Y.death.ALG <- min(1,
                         param$alpha.N.ALG / param$alpha.N.POM,
                         param$alpha.P.ALG / param$alpha.P.POM,
                         param$alpha.C.ALG / param$alpha.C.POM)
param$Y.death.ZOO <- min(1,
                         param$alpha.N.ZOO / param$alpha.N.POM,
                         param$alpha.P.ZOO / param$alpha.P.POM,
                         param$alpha.C.ZOO / param$alpha.C.POM)

# define substances and organisms
NH4  <- c(N      = 1,
          H      = 4 * 1/14,
          charge = 1/14)
NO3  <- c(N      = 1,
          O      = 3 * 16/14,
          charge = -1/14)
N2   <- c(N      = 1)
HPO4 <- c(H      = 1 * 1/31,
          P      = 1,
          O      = 4 * 16/31,
          charge = -2/31)
HCO3 <- c(H      = 1 * 1/12,
          C      = 1,
          O      = 3 * 16/12,
          charge = -1/12)
O2   <- c(O      = 1)
H    <- c(H      = 1,
          charge = 1)
H2O  <- c(H      = 2 * 1,
          O      = 1 * 16)

ALG <- c(O = param$alpha.O.ALG,
         H = param$alpha.H.ALG,
         N = param$alpha.N.ALG,
         P = param$alpha.P.ALG,
         C = param$alpha.C.ALG)
ZOO <- c(O = param$alpha.O.ZOO,
         H = param$alpha.H.ZOO,
         N = param$alpha.N.ZOO,
         P = param$alpha.P.ZOO,
         C = param$alpha.C.ZOO)
POM <- c(O = param$alpha.O.POM,
         H = param$alpha.H.POM,
         N = param$alpha.N.POM,
         P = param$alpha.P.POM,
         C = param$alpha.C.POM)

subst.comp <- list(C.NH4   = NH4,
                   C.NO3   = NO3,
                   C.N2    = N2,
                   C.HPO4  = HPO4,
                   C.HCO3  = HCO3,
                   C.O2    = O2,
                   C.H     = H,
                   C.H2O   = H2O,
                   C.ALG   = ALG,
                   C.ZOO   = ZOO,
                   C.POMD  = POM,
                   C.POMI  = POM,
                   D.POMD  = POM,
                   D.POMI  = POM)


# build composition matrix
alpha <- calc.comp.matrix(subst.comp)
print(round(alpha, 3))


# calculate stoichiometric coefficients
nu.growth.ALG.NH4 <- calc.stoich.coef(alpha       = alpha,
                                      name        = "growth.ALG.NH4",
                                      subst       = c("C.NH4", "C.HPO4", "C.HCO3", "C.O2", "C.H", "C.H2O", "C.ALG"),
                                      subst.norm  = "C.ALG",
                                      nu.norm     = 1.0)
nu.growth.ALG.NO3 <- calc.stoich.coef(alpha       = alpha,
                                      name        = "growth.ALG.NO3",
                                      subst       = c("C.NO3", "C.HPO4", "C.HCO3", "C.O2", "C.H", "C.H2O", "C.ALG"),
                                      subst.norm  = "C.ALG",
                                      nu.norm     = 1.0)                                      
nu.resp.ALG       <- calc.stoich.coef(alpha       = alpha,
                                      name        = "resp.ALG",
                                      subst       = c("C.NH4", "C.HPO4", "C.HCO3", "C.O2", "C.H", "C.H2O", "C.ALG"),
                                      subst.norm  = "C.ALG",
                                      nu.norm     = -1)
nu.death.ALG      <- calc.stoich.coef(alpha       = alpha,
                                      name        = "death.ALG",
                                      subst       = c("C.NH4", "C.HPO4", "C.HCO3", "C.O2", "C.H", "C.H2O", "C.ALG", "C.POMD", "C.POMI"),
                                      subst.norm  = "C.ALG",
                                      nu.norm     = -1,
                                      constraints = list(c("C.ALG"  = param$Y.death.ALG,
                                                           "C.POMD" = 1,
                                                           "C.POMI" = 1),
                                                         c("C.POMD" = -param$f.I,
                                                           "C.POMI" = 1 - param$f.I)))

nu.growth.ZOO <- calc.stoich.coef(alpha       = alpha,
                                  name        = "growth.ZOO",
                                  subst       = c("C.NH4", "C.HPO4", "C.HCO3", "C.O2", "C.H", "C.H2O", "C.ALG", "C.ZOO", "C.POMD", "C.POMI"),
                                  subst.norm  = "C.ZOO",
                                  nu.norm     = 1,
                                  constraints = list(c("C.ALG" = param$Y.ZOO,
                                                       "C.ZOO" = 1),
                                                     c("C.ALG"  = param$f.e,
                                                       "C.POMD" = 1,
                                                       "C.POMI" = 1),
                                                     c("C.POMD" = -param$f.I,
                                                       "C.POMI" = 1 - param$f.I)))
nu.resp.ZOO   <- calc.stoich.coef(alpha       = alpha,
                                  name        = "resp.ZOO",
                                  subst       = c("C.NH4", "C.HPO4", "C.HCO3", "C.O2", "C.H", "C.H2O", "C.ZOO"),
                                  subst.norm  = "C.ZOO",
                                  nu.norm     = -1)                                 
nu.death.ZOO  <- calc.stoich.coef(alpha       = alpha,
                                  name        = "death.ZOO",
                                  subst       = c("C.NH4", "C.HPO4", "C.HCO3", "C.O2", "C.H", "C.H2O", "C.ZOO", "C.POMD", "C.POMI"),
                                  subst.norm  = "C.ZOO",
                                  nu.norm     = -1,
                                  constraints = list(c("C.ZOO"  = param$Y.death.ZOO,
                                                       "C.POMD" = 1,
                                                       "C.POMI" = 1),
                                                     c("C.POMD" = -param$f.I,
                                                       "C.POMI" = 1 - param$f.I)))

nu.nitri <- calc.stoich.coef(alpha      = alpha,
                             name       = "nitri",
                             subst      = c("C.NH4", "C.NO3", "C.O2", "C.H", "C.H2O"),
                             subst.norm = "C.NH4",
                             nu.norm    = -1)

nu.miner.ox.POM       <- calc.stoich.coef(alpha       = alpha,
                                          name        = "miner.ox.POM",
                                          subst       = c("C.NH4", "C.HPO4", "C.HCO3", "C.O2", "C.H", "C.H2O", "C.POMD"),
                                          subst.norm  = "C.POMD",
                                          nu.norm     = -1)
nu.miner.ox.POM.sed   <- calc.stoich.coef(alpha       = alpha,
                                          name        = "miner.ox.POM.sed",
                                          subst       = c("C.NH4", "C.HPO4", "C.HCO3", "C.O2", "C.H", "C.H2O", "D.POMD"),
                                          subst.norm  = "D.POMD",
                                          nu.norm     = -1)
nu.miner.anox.POM.sed <- calc.stoich.coef(alpha       = alpha,
                                          name        = "miner.anox.POM.sed",
                                          subst       = c("C.NH4", "C.NO3", "C.N2", "C.HPO4", "C.HCO3", "C.H", "C.H2O", "D.POMD"),
                                          subst.norm  = "D.POMD",
                                          nu.norm     = -1,
                                          constraints = list(c("C.NO3" = 1,
                                                               "C.N2"  = 1)))

nu.sed.POMD <- calc.stoich.coef(alpha      = alpha,
                                name       = "sed.POMD",
                                subst      = c("C.POMD", "D.POMD"),
                                subst.norm = "C.POMD",
                                nu.norm    = -1)
nu.sed.POMI <- calc.stoich.coef(alpha      = alpha,
                                name       = "sed.POMI",
                                subst      = c("C.POMI", "D.POMI"),
                                subst.norm = "C.POMI",
                                nu.norm    = -1)

nu <- rbind(nu.growth.ALG.NH4,
            nu.growth.ALG.NO3,
            nu.resp.ALG,
            nu.death.ALG,
            nu.growth.ZOO,
            nu.resp.ZOO,
            nu.death.ZOO,
            nu.nitri,
            nu.miner.ox.POM,
            nu.miner.ox.POM.sed,
            nu.miner.anox.POM.sed,
            nu.sed.POMD,
            nu.sed.POMI)
print(round(nu, 3))


# transformation processes
growth.ALG.NH4 <- new(Class  = "process",
                      name   = "Growth of algae with NH4",
                      rate   = expression(k.growth.ALG 
                                           * exp(beta.ALG * (T - T0))
                                           * log((K.I + I0) / (K.I + I0 * exp(-(lambda.1 + lambda.2 * C.ALG) * h.epi)))
                                              / ((lambda.1 + lambda.2 * C.ALG) * h.epi)
                                           * min((C.NH4 + C.NO3) / (K.N + C.NH4 + C.NO3), C.HPO4 / (K.HPO4 + C.HPO4))
                                           * (p.NH4 * C.NH4) / (p.NH4 * C.NH4 + C.NO3)
                                           * C.ALG),
                      stoich = as.list(nu["growth.ALG.NH4", ]))
growth.ALG.NO3 <- new(Class  = "process",
                      name   = "Growth of algae with NO3",
                      rate   = expression(k.growth.ALG 
                                           * exp(beta.ALG * (T - T0))
                                           * log((K.I + I0) / (K.I + I0 * exp(-(lambda.1 + lambda.2 * C.ALG) * h.epi)))
                                              / ((lambda.1 + lambda.2 * C.ALG) * h.epi)
                                           * min((C.NH4 + C.NO3) / (K.N + C.NH4 + C.NO3), C.HPO4 / (K.HPO4 + C.HPO4))
                                           * C.NO3 / (p.NH4 * C.NH4 + C.NO3)
                                           * C.ALG),
                      stoich = as.list(nu["growth.ALG.NO3", ]))
resp.ALG       <- new(Class  = "process",
                      name   = "Respiration of algae",
                      rate   = expression(k.resp.ALG
                                           * exp(beta.ALG * (T - T0))
                                           * C.O2 / (K.O2.resp + C.O2)
                                           * C.ALG),
                      stoich = as.list(nu["resp.ALG", ]))
death.ALG      <- new(Class  = "process",
                      name   = "Death of algae",
                      rate   = expression(k.death.ALG * C.ALG),
                      stoich = as.list(nu["death.ALG", ]))

growth.ZOO <- new(Class  = "process",
                  name   = "Growth of zooplankton",
                  rate   = expression(k.growth.ZOO
                                       * exp(beta.ZOO * (T - T0))
                                       * C.O2 / (K.O2.ZOO + C.O2)
                                       * C.ALG
                                       * C.ZOO),
                  stoich = as.list(nu["growth.ZOO", ]))
resp.ZOO   <- new(Class  = "process",
                  name   = "Respiration of zooplankton",
                  rate   = expression(k.resp.ZOO
                                       * exp(beta.ZOO * (T - T0))
                                       * C.O2 / (K.O2.resp + C.O2)
                                       * C.ZOO),
                  stoich = as.list(nu["resp.ZOO", ]))
death.ZOO  <- new(Class  = "process",
                  name   = "Death of zooplankton",
                  rate   = expression(k.death.ZOO * C.ZOO),
                  stoich = as.list(nu["death.ZOO", ]))

nitri <- new(Class  = "process",
             name   = "Nitrification",
             rate   = expression(k.nitri
                                  * exp(beta.BAC * (T - T0))
                                  * min(C.NH4 / (K.NH4.nitri + C.NH4), C.O2 / (K.O2.nitri + C.O2))),
             stoich = as.list(nu["nitri", ]))

miner.ox.POM       <- new(Class  = "process",
                          name   = "Oxic mineralization",
                          rate   = expression(k.miner.ox.POM
                                               * exp(beta.BAC * (T - T0))
                                               * C.O2 / (K.O2.miner + C.O2)
                                               * C.POMD),
                          stoich = as.list(nu["miner.ox.POM", ]))
miner.ox.POM.sed   <- new(Class  = "process",
                          name   = "Oxic mineralization in sediment",
                          rate   = expression(k.miner.ox.POM.sed
                                               * exp(beta.BAC * (T - T0))
                                               * C.O2 / (K.O2.miner + C.O2)
                                               * D.POMD / (K.POM.miner.sed + D.POMD)),
                          stoich = as.list(nu["miner.ox.POM.sed", ]),
                          pervol = F)
miner.anox.POM.sed <- new(Class  = "process",
                          name   = "Anoxic mineralization in sediment",
                          rate   = expression(k.miner.anox.POM.sed
                                               * exp(beta.BAC * (T - T0))
                                               * C.NO3 / (K.NO3.miner + C.NO3)
                                               * (D.POMD / (K.POM.miner.sed + D.POMD))^2),
                          stoich = as.list(nu["miner.anox.POM.sed", ]),
                          pervol = F)

sed.POMD <- new(Class  = "process",
                name   = "Sedimentation of Dead Particulate Organic Matter",
                rate   = expression(v.sed.POM / h.hypo * C.POMD),
                stoich = as.list(nu["sed.POMD", ]))
sed.POMI <- new(Class  = "process",
                name   = "Sedimentation of Inert Particulate Organic Matter",
                rate   = expression(v.sed.POM / h.hypo * C.POMI),
                stoich = as.list(nu["sed.POMI", ]))


# environmental conditions
cond.epi  <- list(I0       = expression(0.5 * (I0.min + I0.max) 
                                         + 0.5 * (I0.max - I0.min)
                                          * cos(2*pi / 365.25 * (t - t.max))),
                  T        = expression(0.5 * (T.min + T.max) 
                                         + 0.5 * (T.max - T.min)
                                          * cos(2*pi / 365.25 * (t - t.max))),
                  C.O2.sat = expression(exp(7.7117 - 1.31403 * log(T + 45.93))
                                         * p / 101325))
cond.hypo <- list(I0       = 0,
                  T        = 5)
cond.gen  <- list(Kz       = expression(0.5 * (Kz.summer + Kz.winter)
                                         - 0.5 * (Kz.winter - Kz.summer)
                                          * sign(cos(2*pi / 365.25 * (t - t.max)) + 0.4)))

t        <- 1:365
Kz       <- numeric(0)
I0       <- numeric(0)
T        <- numeric(0)
C.O2.sat <- numeric(0)
for(i in 1:length(t)) {
  Kz[i]       <- eval(cond.gen$Kz,       envir= c(param, t= t[i]))
  I0[i]       <- eval(cond.epi$I0,       envir= c(param, t= t[i]))
  T[i]        <- eval(cond.epi$T,        envir= c(param, t= t[i]))
  C.O2.sat[i] <- eval(cond.epi$C.O2.sat, envir= c(param, T= T[i]))
}

par.def <- par(no.readonly= TRUE)
par(mfrow = c(2,2),
    xaxs  = "i",
    yaxs  = "i", 
    mar   = c(4.5, 4.5, 2, 1.5) + 0.1)
plot(t, Kz,       ylim= c(0,25))
plot(t, I0,       type= "l")
plot(t, T,        type= "l")
plot(t, C.O2.sat, type= "l")
par(par.def)


# reactors
epilimnion  <- new(Class            = "reactor",
                   name             = "Epi",
                   volume.ini       = expression(A * h.epi),
                   conc.pervol.ini  = list(C.NH4  = expression(C.NH4.ini),
                                           C.NO3  = expression(C.NO3.ini),
                                           C.HPO4 = expression(C.HPO4.ini),
                                           C.O2   = expression(C.O2.ini),
                                           C.ALG  = expression(C.ALG.ini),
                                           C.ZOO  = expression(C.ZOO.ini),
                                           C.POMD = expression(C.POMD.ini),
                                           C.POMI = expression(C.POMI.ini)),
                   input            = list(C.O2   = expression(v.ex.O2 * A
                                                                * (C.O2.sat - C.O2))),
                   inflow           = expression(Q.in * 86400),
                   inflow.conc      = list(C.NH4  = 0.0,
                                           C.NO3  = expression(C.NO3.in),
                                           C.HPO4 = expression(C.HPO4.in),
                                           C.O2   = expression(C.O2.in),
                                           C.ALG  = 0.0,
                                           C.ZOO  = 0.0,
                                           C.POMD = 0.0,
                                           C.POMI = 0.0),
                   outflow          = expression(Q.in * 86400),
                   cond             = cond.epi,
                   processes        = list(growth.ALG.NH4, growth.ALG.NO3, resp.ALG, death.ALG,
                                           growth.ZOO, resp.ZOO, death.ZOO,
                                           nitri, miner.ox.POM))
hypolimnion <- new(Class            = "reactor",
                   name             = "Hypo",
                   volume.ini       = expression(A * h.hypo),
                   conc.pervol.ini  = list(C.NH4  = expression(C.NH4.ini),
                                           C.NO3  = expression(C.NO3.ini),
                                           C.HPO4 = expression(C.HPO4.ini),
                                           C.O2   = expression(C.O2.ini),
                                           C.ALG  = expression(C.ALG.ini),
                                           C.ZOO  = expression(C.ZOO.ini),
                                           C.POMD = expression(C.POMD.ini),
                                           C.POMI = expression(C.POMI.ini)),
                   area             = expression(A),
                   conc.perarea.ini = list(D.POMD = expression(D.POMD.ini),
                                           D.POMI = expression(D.POMI.ini)),
                   cond             = cond.hypo,
                   processes        = list(resp.ALG, death.ALG,
                                           growth.ZOO, resp.ZOO, death.ZOO,
                                           nitri,
                                           miner.ox.POM, miner.ox.POM.sed, miner.anox.POM.sed,
                                           sed.POMD, sed.POMI))

# link
metalimnion <- new(Class     = "link",
                   name      = "Meta",
                   from      = "Epi",
                   to        = "Hypo",
                   qadv.spec = list(C.POMD = expression(v.sed.POM * A),
                                    C.POMI = expression(v.sed.POM * A)),
                   qdiff.gen = expression((A / h.meta) * Kz))


# system
lake <- new(Class    = "system",
            name     = "Lake",
            reactors = list(epilimnion, hypolimnion),
            links    = list(metalimnion),
            cond     = cond.gen,
            param    = param,
            t.out    = seq(0, 730, by= 1))


# run simulation and plot result
ptm <- proc.time()
res <- calcres(lake)
print(proc.time() - ptm)
ptm <- proc.time()
res.rk4 <- calcres(lake, method= "rk4")
print(proc.time() - ptm)

plotres(res      = res,
        colnames = list(c("C.NH4.Epi",   "C.NH4.Hypo"),
                        c("C.NO3.Epi",   "C.NO3.Hypo"),
                        c("C.HPO4.Epi",  "C.HPO4.Hypo"),
                        c("C.O2.Epi",    "C.O2.Hypo"),
                        c("C.ALG.Epi",   "C.ALG.Hypo"),
                        c("C.ZOO.Epi",   "C.ZOO.Hypo"),
                        c("C.POMD.Epi",  "C.POMD.Hypo"),
                        c("C.POMI.Epi",  "C.POMI.Hypo"),
                        c("D.POMD.Hypo", "D.POMI.Hypo")))
plotres(res      = list(lsoda= res, rk4= res.rk4),  # compare numerical solutions
        colnames = list(c("C.NH4.Epi",   "C.NH4.Hypo"),
                        c("C.NO3.Epi",   "C.NO3.Hypo"),
                        c("C.HPO4.Epi",  "C.HPO4.Hypo"),
                        c("C.O2.Epi",    "C.O2.Hypo"),
                        c("C.ALG.Epi",   "C.ALG.Hypo"),
                        c("C.ZOO.Epi",   "C.ZOO.Hypo"),
                        c("C.POMD.Epi",  "C.POMD.Hypo"),
                        c("C.POMI.Epi",  "C.POMI.Hypo"),
                        c("D.POMD.Hypo", "D.POMI.Hypo")))