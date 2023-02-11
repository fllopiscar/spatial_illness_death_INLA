library(INLA)
library(rgdal)
library(rgeos)
library(dplyr)
library(rlang)
library(stringr)
library(spdep)

#Adjacency matrix
load("map_depart.R")
Valencia.nb <- poly2nb(mapa_depart)
adj <- Valencia.nb
W <- as(nb2mat(adj, style = "B"), "sparseMatrix")

#Data preparation
joint.data <- data_preparation()

#Define latent effects model
source("multiLeroux3d function.R")
multiLeroux3d.model <- inla.rgeneric.define(inla.rgeneric.multiLeroux3d.model, W = W, v=7,d.elem=1)
#Formula 
formula.model.re.leroux=Y~-1+beta0+sex1+sex2+age1+age2+sex23+age23+
                             f(dep.aux, model=multiLeroux3d.model)
#INLA call
multistate.multiLeroux.wishart1.full= inla(formula.model.re.leroux, family=c("weibullsurv", "weibullsurv", "weibullsurv"), 
                                    data=joint.data, control.compute = list(dic=TRUE, config=TRUE),
                                    control.predictor = list(compute = TRUE), verbose=T)

#Summary and plot of random effects
summary(multistate.multiLeroux.wishart1.full)

mapa_depart$multiLeroux_whishart1_1 <- multistate.multiLeroux.wishart1.full$summary.random$h.area[1:24, "mean"]
mapa_depart$multiLeroux_whishart1_2 <- multistate.multiLeroux.wishart1.full$summary.random$h.area[25:48, "mean"]
mapa_depart$multiLeroux_whishart1_3 <- multistate.multiLeroux.wishart1.full$summary.random$h.area[49:72, "mean"]

spplot(mapa_depart, c("multiLeroux_whishart1_1","multiLeroux_whishart1_2","multiLeroux_whishart1_3"),
       col.regions=rev(Paleta(16)))
