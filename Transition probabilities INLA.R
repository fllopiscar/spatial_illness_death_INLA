#We draw from the joint posterior distribution
for(i in 1:20){
post_sample.1  <- inla.posterior.sample(n=100, multistate.multiLeroux.wishart1.full, seed = 12345+i)

post_samplef.1 <- inla.posterior.sample.eval(function(...) c(beta01, beta02, beta03, sex1, sex2, age1, age2, sex23, age23,
                                                             theta[1:3], dep.aux), post_sample.1)
rm(post_sample.1)

sample.1 <- t(post_samplef.1)
colnames(sample.1) <- c("beta12", "beta13","beta23", "sex12", "sex13", "age12", "age13", "sex23", "age23", 
                                 "alpha12", "alpha13", "alpha23", paste("dep", 1:72, sep=""))
sample.1 <- data.frame(sample.1)
rm(post_samplef.1)
if(i==1){sample.posterior <- sample.1}
else{sample.posterior <- rbind(sample.posterior, sample.1)}
}

#Simulated sample
sample.posterior.wishart1 <- sample.posterior 

############################################################
library(INLA)
edad_i <- datos_cohorte$edad_i
attach(sample.posterior.wishart1)

############################################
# F12 - Cumulative incidence of refracture #
############################################

cif1.inla <- function(t,x){
  media<- numeric(length(t))
  quantiles<- matrix(nrow=2,ncol=length(t))
  int<-   numeric(length(alpha12))
  
  x[2]<- x[2]-mean(edad_i)
  
  effect.dep.12 <- eval(parse(text=paste("dep", x[3], sep="")))
  effect.dep.13 <- eval(parse(text=paste("dep", 24+x[3], sep="")))
  effect.dep.23 <- eval(parse(text=paste("dep", 2*24+x[3], sep="")))
  
  
  exp12<-exp(x[1]*sex12+x[2]*age12+effect.dep.12)
  exp13<-exp(x[1]*sex13+x[2]*age13+effect.dep.13)
  lambda12 <- exp(beta12)
  lambda13 <- exp(beta13)
  
  ext<- lambda12*alpha12*exp12
  
  for(i in 1:length(t)){
    
    for(j in 1:length(alpha12)){
      
      ff<-function(u){u^(alpha12[j]-1)*exp(-lambda12[j]*u^alpha12[j]*exp12[j]-lambda13[j]*u^alpha13[j]*exp13[j])}
      int[j]<-integrate(ff,0,t[i],rel.tol = 1e-8)[[1]]
    }
    
    media[i]<- mean(ext*int)
    quantiles[,i] <- quantile(ext*int, probs=c(0.025,0.975))
  }
  
  list(media,quantiles)
  
}


################################
# P11 - Event-free probability #
################################

p11.func <- function(s,t,x){
  media<- numeric(length(t))
  quantiles<- matrix(nrow=2,ncol=length(t))
  x[2]<- x[2]-mean(edad_i)
  int<-   numeric(length(alpha12))
  
  effect.dep.12 <- eval(parse(text=paste("dep", x[3], sep="")))
  effect.dep.13 <- eval(parse(text=paste("dep", 24+x[3], sep="")))
  
  exp12<-exp(x[1]*sex12+x[2]*age12+effect.dep.12)
  exp13<-exp(x[1]*sex13+x[2]*age13+effect.dep.13)
  lambda12 <- exp(beta12)
  lambda13 <- exp(beta13)
  
  for(i in 1:length(t)){
    aux<- exp(-lambda12*(t[i]^alpha12-s^alpha12)*exp12-lambda13*(t[i]^alpha13-s^alpha13)*exp13)
    media[i]<-mean(aux)
    quantiles[,i] <- quantile(aux, probs=c(0.025,0.975))
  }
  list(media, quantiles)
}


#################################
# P13 - Total death probability #
#################################

p13.func<- function(s,t,x){
  
  media<- numeric(length(t))
  quantiles<- matrix(nrow=2,ncol=length(t))
  x[2]<- x[2]-mean(edad_i)
  int<-   numeric(length(alpha12))
  
  effect.dep.12 <- eval(parse(text=paste("dep", x[3], sep="")))
  effect.dep.13 <- eval(parse(text=paste("dep", 24+x[3], sep="")))
  effect.dep.23 <- eval(parse(text=paste("dep", 2*24+x[3], sep="")))
  
  exp12<-exp(x[1]*sex12+x[2]*age12+effect.dep.12)
  exp13<-exp(x[1]*sex13+x[2]*age13+effect.dep.13)
  exp23<-exp(x[1]*sex23+x[2]*age23+effect.dep.23)
  lambda12 <- exp(beta12)
  lambda13 <- exp(beta13)
  lambda23 <- exp(beta23)
  
  ext<- lambda12*alpha12*exp12*exp(lambda12*s^alpha12*exp12+lambda13*s^alpha13*exp13)
  
  for(i in 1:length(t)){
    p11<- exp(-lambda12*(t[i]^alpha12-s^alpha12)*exp12-lambda13*(t[i]^alpha13-s^alpha13)*exp13)
    
    for(j in 1:length(alpha12)){
      
      ff<-function(u){u^(alpha12[j]-1)*exp(-lambda12[j]*u^alpha12[j]*exp12[j]-lambda13[j]*u^alpha13[j]*exp13[j]-
                                             lambda23[j]*(t[i]-u)^alpha23[j]*exp23[j])}
      int[j]<-integrate(ff,s,t[i],rel.tol = 1e-8)[[1]]
    }
    aux<- 1-p11-int*ext
    media[i]<- mean(aux)
    quantiles[,i] <- quantile(aux, probs=c(0.025,0.975))
    
  }
  
  list(media,quantiles)
  
}

############################################
# P12 - Refractured and alive probability  #
############################################

p12.func<- function(s,t,x){
  
  media<- numeric(length(t))
  quantiles<- matrix(nrow=2,ncol=length(t))
  x[2]<- x[2]-mean(edad_i)
  int<-   numeric(length(alpha12))
  
  effect.dep.12 <- eval(parse(text=paste("dep", x[3], sep="")))
  effect.dep.13 <- eval(parse(text=paste("dep", 24+x[3], sep="")))
  effect.dep.23 <- eval(parse(text=paste("dep", 2*24+x[3], sep="")))
  
  exp12<-exp(x[1]*sex12+x[2]*age12+effect.dep.12)
  exp13<-exp(x[1]*sex13+x[2]*age13+effect.dep.13)
  exp23<-exp(x[1]*sex23+x[2]*age23+effect.dep.23)
  lambda12 <- exp(beta12)
  lambda13 <- exp(beta13)
  lambda23 <- exp(beta23)
  
  ext<- lambda12*alpha12*exp12*exp(lambda12*s^alpha12*exp12+lambda13*s^alpha13*exp13)
  
  for(i in 1:length(t)){
    
    for(j in 1:length(alpha12)){
      
      ff<-function(u){u^(alpha12[j]-1)*exp(-lambda12[j]*u^alpha12[j]*exp12[j]-lambda13[j]*u^alpha13[j]*exp13[j]-
                                             lambda23[j]*(t[i]-u)^alpha23[j]*exp23[j])}
      int[j]<-integrate(ff,s,t[i],rel.tol = 1e-8)[[1]]
    }
    
    media[i]<- mean(ext*int)
    quantiles[,i] <- quantile(ext*int, probs=c(0.025,0.975))
  }
  
  list(media, quantiles)
  
}




############################################
# P23 - Death after refracture probability #
############################################

p23.func <- function(s,t,t12,x){
  media<- numeric(length(t))
  quantiles<- matrix(nrow=2,ncol=length(t))
  x[2]<- x[2]-mean(edad_i)
  int<-   numeric(length(alpha12))
  
  effect.dep.23 <- eval(parse(text=paste("dep", 2*24+x[3], sep="")))
  exp23<-exp(x[1]*sex23+x[2]*age23+effect.dep.23)
  lambda23 <- exp(beta23)
  
  for(i in 1:length(t)){
    aux<- 1-exp(-lambda23*((t[i]-t12)^alpha23-(s-t12)^alpha23)*exp23)
    media[i]<-mean(aux)
    quantiles[,i] <- quantile(aux, probs=c(0.025,0.975))
  }
  list(media,quantiles)
  
}
