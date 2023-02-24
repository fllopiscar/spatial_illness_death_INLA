library(INLA)
library(rgdal)
library(rgeos)
library(dplyr)
library(rlang)
library(stringr)
library(spdep)

library(rjags)
library(snow)
library(mcmcplots)

load("mapa_depart.R")
Valencia.nb <- poly2nb(mapa_depart)
adj <- Valencia.nb

source("Functions/data_preparation.R")

mapa_full <- mapa_depart
#list_dep <- c(20,21,22)
list_dep <- c(12,13,15,16)
list_dep <- c(6,8,9,14,23)
#mapa_depart <- mapa_depart[mapa_depart$dep%in% c(12,13,15,16),]
mapa_depart <- mapa_full[mapa_full$dep%in% list_dep,]
plot(mapa_depart)

################# INLA ####################

joint.data <- data_preparation(N=5000, all=T, list_dep=list_dep)

formula1 = Y~-1+beta0+sex1+sex2+age1+age2+sex23+age23+f(dep.aux, model="iid3d", n=3*length(list_dep),
                                                        constr = F)

model1 = inla(formula1, family=c("weibullsurv", "weibullsurv", "weibullsurv"), data=joint.data,
                   control.predictor = list(compute = TRUE),
                   verbose=F)

model1$summary.random
summary(model1)

################# JAGS ######################

source("Functions/JAGS model iid3d.R")
source("Functions/data_preparation_jags.R")

d.jags <- data_preparation_jags(N=10000, all=T, list_dep=list_dep)

p.jags <- c("beta", "alpha", "beta0", "a", "rho", "Sigma")

coda.samples.wrapper <- function(j)
{ 
  m1 = jags.model(data=d.jags, file=textConnection(weib.ms.iid3d.model), 
                  inits=list(.RNG.name="base::Wichmann-Hill", .RNG.seed=j), n.chains=1)
  update(m1, 4000)
  coda.samples(m1, variable.names=p.jags, n.iter=40000, thin=1) 
}

n.chains=3
snow.start.time = proc.time()
cl <- makeCluster(n.chains, "SOCK")
##Make sure the rjags library is loaded in each worker
clusterEvalQ(cl, library(rjags))
##Send data to workers, then fit models. One disadvantage of this
##parallelization is that you lose the ability to watch the progress bar.
clusterExport(cl, list("d.jags", "weib.ms.iid3d.model", "p.jags"))
samples = clusterApply(cl, 1:3, coda.samples.wrapper)
##Reorganize 'par.samples' so that it is recognizeable as an 'mcmc.list' object
for(i in 1:length(samples)) {samples[[i]] <- samples[[i]][[1]] }
class(samples) <- "mcmc.list"
stopCluster(cl)
snow.end.time = proc.time()
snow.dtime = snow.end.time - snow.start.time
snow.dtime

mcmcplot(samples)
beep()


sum_samples <- summary(samples)
mapa_depart$jags_12 <- sum_samples$statistics[10:13,1]
mapa_depart$jags_13 <- sum_samples$statistics[14:17,1]
mapa_depart$jags_23 <- sum_samples$statistics[18:21,1]

mapa_depart$inla_12 <- model1$summary.random$dep.aux[1:length(list_dep), 2]
mapa_depart$inla_13 <- model1$summary.random$dep.aux[length(list_dep)+1:length(list_dep), 2]
mapa_depart$inla_23 <- model1$summary.random$dep.aux[2*length(list_dep)+1:length(list_dep), 2]

library(viridisLite)
library(RColorBrewer)

Paleta2<-colorRampPalette(c("Darkblue","White","Red3"))(16)

spplot(mapa_depart,c("inla_12", "inla_13","inla_23"), col.regions=rev(rainbow(16, start=0, end=0.3)))

Paleta2<-colorRampPalette(c("#20613c","#429667","lightyellow","Red3","Red4"))(16)
spplot(mapa_depart,c("inla_12", "inla_13","inla_23"), col.regions=Paleta2)


mapa_full$inla_12 <- model_full$summary.random$dep.aux$mean[1:24]
mapa_full$inla_13 <- model_full$summary.random$dep.aux$mean[24+1:24]
mapa_full$inla_23 <- model_full$summary.random$dep.aux$mean[2*24+1:24]

#Paleta terreno (similar)
Paleta2<-colorRampPalette(c("#20613c","#429667","lightyellow","tan3","#560909"))(16)
spplot(mapa_full, c("inla_12", "inla_13", "inla_23") , col.regions=Paleta2)
Paleta2<-colorRampPalette(c("#20613c","#52a878","lightyellow","tan3","#560909"))(16)
#Paleta terreno (similar)
Paleta2<-colorRampPalette(c("Blue4","lightyellow","Red3"))(16)
spplot(mapa_full, c("inla_12", "inla_13", "inla_23") , col.regions=Paleta2)

#spplot(mapa_full, c("inla_12", "inla_13", "inla_23") , col.regions=rev(hcl.colors(16, "Reds")))
spplot(mapa_full, c("inla_12", "inla_13", "inla_23") , col.regions=rev(heat.colors(16)))


spplot(mapa_depart,c("jags_12", "inla_12"), col.regions=Paleta2)
spplot(mapa_depart,c("jags_13", "inla_13"), col.regions=Paleta2)
spplot(mapa_depart,c("jags_23", "inla_23"), col.regions=Paleta2)

spplot(mapa_depart, c("inla_12", "jags_12", "inla_13", "jags_13", "inla_23","jags_23") , col.regions=Paleta2)

gelman.diag(samples, multivariate = FALSE)

mapa_depart_aux <- mapa_depart
mapa_depart_aux$inla_vs_jags_12 <- mapa_depart$inla_12-mapa_depart$jags_12
mapa_depart_aux$inla_vs_jags_13 <- mapa_depart$inla_13-mapa_depart$jags_13
mapa_depart_aux$inla_vs_jags_23 <- mapa_depart$inla_23-mapa_depart$jags_23

spplot(mapa_depart_aux, c("inla_vs_jags_12", "inla_vs_jags_13", "inla_vs_jags_23") , col.regions=Paleta2)


summary(abs(c(mapa_depart_aux$inla_vs_jags_12,mapa_depart_aux$inla_vs_jags_13,mapa_depart_aux$inla_vs_jags_23)))

