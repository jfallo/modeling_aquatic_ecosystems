library("ecosim")
library("stoichcalc")
library("deSolve")


# define model parameters
param <- list(k.gro.ALG   = 0.5,     # 1/d
              k.death.ALG = 0.1,     # 1/d
              K.HPO4      = 0.002,   # gP/m3
              alpha.P.ALG = 0.003,   # gP/gDM
              A           = 5e+006,  # m2
              h.epi       = 5,       # m
              Q.in        = 5,       # m3/s
              C.HPO4.in   = 0.04,    # gP/m3
              C.HPO4.ini  = 0.004,   # gP/m3
              C.ALG.ini   = 0.1)     # gDM/m3

# define transformation processes
gro.ALG <- new(Class  = "process",
               name   = "Growth of algae",
               rate   = expression(k.gro.ALG *C.HPO4/(K.HPO4 + C.HPO4) * C.ALG),
               stoich = list(C.ALG  = expression(1),              # gDM/gDM
                             C.HPO4 = expression(-alpha.P.ALG)))  # gP/gDM
death.ALG <- new(Class  = "process",
                 name   = "Death of algae",
                 rate   = expression(k.death.ALG * C.ALG),
                 stoich = list(C.ALG = expression(-1)))  # gDM/gDM

# define reactor to describe the epilimnion of the lake
epilimnion <- new(Class           = "reactor",
                  name            = "Epilimnion",
                  volume.ini      = expression(A * h.epi),
                  conc.pervol.ini = list(C.HPO4 = expression(C.HPO4.ini),  # gP/m3
                                         C.ALG = expression(C.ALG.ini)),   # gDM/m3
                  inflow          = expression(Q.in * 86400),                # m3/d
                  inflow.conc     = list(C.HPO4 = expression(C.HPO4.in),
                                         C.ALG = 0),
                  outflow         = expression(Q.in * 86400),
                  processes       = list(gro.ALG, death.ALG))

# define system consisting of a single reactor
system.lake_phytoplankton <- new(Class    = "system",
                                 name     = "Lake",
                                 reactors = list(epilimnion),
                                 param    = param,
                                 t.out    = seq(0, 365, by= 1))

# simulate and plot results
res <- calcres(system.lake_phytoplankton)
plotres(res, colnames= list("C.HPO4", "C.ALG"))

# change to 4 year analysis
system.lake_phytoplankton@t.out <- seq(0, 4*365, by= 1)
res.4y <- calcres(system.lake_phytoplankton)
plotres(res.4y, colnames= list("C.HPO4", "C.ALG"))


# extend model for seasonally varying driving forces
system.lake_phytoplankton.ext <- system.lake_phytoplankton

# update algae growth
gro.ALG.ext <- new(Class = "process",
                   name  = "Growth of algae extended",
                   rate  = expression(k.gro.ALG 
                                       * exp(beta.ALG * (T - T0))
                                       * C.HPO4 / (K.HPO4 + C.HPO4)
                                       * log((K.I + I0) / (K.I + I0*exp(-(lambda.1 + lambda.2*C.ALG) * h.epi)))
                                          / ((lambda.1 + lambda.2*C.ALG) * h.epi)
                                       * C.ALG),
                   stoich  = list(C.ALG  = 1,                          # gDM/gDM
                                  C.HPO4 = expression(-alpha.P.ALG)))  # gP/gDM
epilimnion@processes <- list(gro.ALG.ext, death.ALG)

# update epilimnion for time-dependent environmental conditions
epilimnion@cond <- list(I0 = expression(0.5*(I0.min + I0.max)
                                         + 0.5*(I0.max - I0.min) * cos(2*pi/365.25 * (t - t.max))),  # W/m2
                        T  = expression(0.5*(T.min + T.max)
                                         + 0.5*(T.max - T.min) * cos(2*pi/365.25 * (t - t.max))))    # degC
system.lake_phytoplankton.ext@reactors <- list(epilimnion)

# extend model parameters
param <- c(param,
           list(beta.ALG = 0.046,  # 1/degC
                T0       = 20,     # degC
                K.I      = 30,     # W/m2
                lambda.1 = 0.10,   # 1/m
                lambda.2 = 0.10,   # m2/gDM
                t.max    = 230,    # d
                I0.min   = 25,     # W/m2
                I0.max   = 225,    # W/m2
                T.min    = 5,      # degC
                T.max    = 25))    # degC
param$k.gro.ALG <- 0.8
system.lake_phytoplankton.ext@param <- param

# simulate and plot results
res.ext <- calcres(system.lake_phytoplankton.ext)
plotres(res.ext, colnames= list("C.HPO4", "C.ALG"))
plotres(res.ext[1:365,], colnames= list("C.HPO4", "C.ALG"))


# sensitivity analysis
sens.res <- calcsens(system          = system.lake_phytoplankton.ext,
                     param.sens      = c("C.HPO4.in", "k.gro.ALG", "k.death.ALG", "K.HPO4"),
                     scaling.factors = c(1, 0.5, 2))
plotres(sens.res)