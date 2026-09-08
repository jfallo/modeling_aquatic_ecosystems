library("ecosim")

# set parameters
param <- list(k.gro.ALG   = 0.5,     # 1/d
              k.gro.ZOO   = 0.4,     # gDM/m3d
              k.death.ALG = 0.1,     # 1/d
              k.death.ZOO = 0.05,    # 1/d
              K.HPO4      = 0.002,   # gP/m3
              Y.ZOO       = 0.2,     # gDM/gDM
              alpha.P.ALG = 0.003,   # gP/gDM
              A           = 5e+006,  # m2
              h.epi       = 5,       # m
              Q.in        = 5,       # m3/s
              C.HPO4.in   = 0.04,    # gP/m3
              C.HPO4.ini  = 0.04,    # gP/m3
              C.ALG.ini   = 0.1,     # gDM/m3
              C.ZOO.ini   = 0.1)     # gDM/m3

# transformation processes
gro.ALG <- new(Class  = "process",
               name   = "Growth of algae",
               rate   = expression(k.gro.ALG * C.HPO4/(K.HPO4 + C.HPO4) * C.ALG),
               stoich = list(C.ALG  = expression(1),
                             C.HPO4 = expression(-alpha.P.ALG)))
death.ALG <- new(Class  = "process",
                 name   = "Death of algae",
                 rate   = expression(k.death.ALG * C.ALG),
                 stoich = list(C.ALG = expression(-1)))

gro.ZOO <- new(Class  = "process",
               name   = "Growth of zooplankton",
               rate   = expression(k.gro.ZOO * C.ALG * C.ZOO),
               stoich = list(C.ALG = expression(-1/Y.ZOO),
                             C.ZOO = expression(1)))
death.ZOO <- new(Class  = "process",
                 name   = "Death of zooplankton",
                 rate   = expression(k.death.ZOO * C.ZOO),
                 stoich = list(C.ZOO = expression(-1)))

# reactor described by lake epilimnion
epilimnion <- new(Class           = "reactor",
                  name            = "Epilimnion",
                  volume.ini      = expression(A * h.epi),
                  conc.pervol.ini = list(C.HPO4 = expression(C.HPO4.ini),
                                         C.ALG  = expression(C.ALG.ini),
                                         C.ZOO  = expression(C.ZOO.ini)),
                  inflow          = expression(Q.in * 86400),
                  inflow.conc     = list(C.HPO4 = expression(C.HPO4.in),
                                         C.ALG  = 0,
                                         C.ZOO  = 0),
                  outflow         = expression(Q.in * 86400),
                  processes       = list(gro.ALG, death.ALG, gro.ZOO, death.ZOO))

# single-reactor system
system.lake_phytozooplankton <- new(Class    = "system",
                                    name     = "Lake",
                                    reactors = list(epilimnion),
                                    param    = param,
                                    t.out    = seq(0, 4*365, by= 1))

# simulate and plot results
res <- calcres(system.lake_phytozooplankton)
plotres(res, colnames= list("C.ALG", "C.ZOO", "C.HPO4"))
plotres(res, colnames= c("C.ALG", "C.ZOO", "C.HPO4"))


# extend model for seasonally varying driving forces
system.lake_phytozooplankton.ext <- system.lake_phytozooplankton

# extend model parameters
param <- c(param,
           list(beta.ALG = 0.046,  # 1/degC
                beta.ZOO = 0.08,   # 1/degC
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
system.lake_phytozooplankton.ext@param <- param

# extend growth processes
gro.ALG.ext <- new(Class  = "process",
                   name   = "Growth of algae",
                   rate   = expression(k.gro.ALG 
                                        * exp(beta.ALG * (T - T0))
                                        * C.HPO4/(K.HPO4 + C.HPO4) 
                                        * log((K.I + I0) / (K.I + I0*exp(-(lambda.1 + lambda.2*C.ALG) * h.epi)))
                                           / ((lambda.1 + lambda.2*C.ALG) * h.epi)
                                        * C.ALG),
                   stoich = list(C.ALG  = expression(1),
                                 C.HPO4 = expression(-alpha.P.ALG)))
gro.ZOO.ext <- new(Class  = "process",
                   name   = "Growth of zooplankton",
                   rate   = expression(k.gro.ZOO 
                                        * exp(beta.ZOO * (T - T0))
                                        * C.ALG 
                                        * C.ZOO),
                   stoich = list(C.ALG = expression(-1/Y.ZOO),
                                 C.ZOO = expression(1)))

# update reactor
epilimnion@processes <- list(gro.ALG.ext, death.ALG, gro.ZOO.ext, death.ZOO)
epilimnion@cond <- list(I0 = expression(0.5*(I0.min + I0.max)
                                         + 0.5*(I0.max - I0.min) * cos(2*pi/365.25 * (t - t.max))),
                        T  = expression(0.5*(T.min + T.max)
                                         + 0.5*(T.max - T.min) * cos(2*pi/365.25 * (t - t.max))))

# update lake system
system.lake_phytozooplankton.ext@reactors <- list(epilimnion)

# plot environmental conditions
t <- seq(1, 4*365)
I0 <- numeric(0)
T <- numeric(0)
for(i in 1:length(t)) {
    I0[i] <- eval(epilimnion@cond$I0, envir= c(param, t= t[i]))
    T[i] <- eval(epilimnion@cond$T, envir= c(param, t= t[i]))
}
par(mfrow= c(1,2), xaxs= "i", yaxs= "i", mar= c(4.5, 4.5, 2, 1.5) + 0.1)
plot(t, I0, type= "l")
plot(t, T, type= "l")

# simulate and plot results
res.ext <- calcres(system.lake_phytozooplankton.ext)
plotres(res.ext, colnames= list("C.HPO4", c("C.ALG", "C.ZOO")))


# compare models
plotres(res = list(const= res, dyn= res.ext),
        colnames = list("C.HPO4", "C.ALG", "C.ZOO"))
