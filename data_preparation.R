data_preparation <- function(){

load("data_cohort.RData")
N <- length(data_cohort)

attach(data_cohort)

#Construction of time-to-event and censoring matrixes

ini <- findex  #discharge after index fracture
end  <- rep(as.Date("2016-12-31", format = "%Y-%m-%d"),length(inic))   #end of study

event <- rep(0,length(findex))                        #no event
event[is.na(fdef)==FALSE]<- 2                         #death
event[is.na(frefrac)==FALSE]<- 1                      #refracture
event[is.na(fdef)==FALSE & is.na(frefrac)==FALSE]<- 3 #death after refracture

c <- matrix(c(as.integer(event== 1), as.integer(event == 2), as.integer(event == 3)), ncol=3)
c[,1][c[,3]==1]<-1

#Time: first column, time from initial to refracture, death or end of study.
#Second column time from initial to death after refracture or end of study..

time <- matrix(0,ncol=2, nrow=length(fdef)) 
time[,1][event==0]           <- (end-findex)[event==0]   
time[,1][event==1| event==3] <- (frefrac-findex)[event==1 | event==3]
time[,1][event==2]           <- (fdef-findex)[event==2]

time[,2][event!=3] <- (end-findex)[event!=3]
time[,2][event==3] <- (fdef-findex)[event==3]

time<-time/365.25


n <- nrow(time)
refrac<- sum(c[,1]==1)

#Times to refracture and to death without refracture (F->R and F->D)
timeC1 <- inla.surv(time=c(time[,1], rep(NA,n), rep(NA,refrac)), event=c(c[,1], rep(NA,n), rep(NA,refrac)))
timeC2 <- inla.surv(time=c(rep(NA,n), time[,1], rep(NA,refrac)), event=c(rep(NA,n), c[,2], rep(NA,refrac)))

#Time after refracture (R->D, only for those refractured) 
time.aux<- (time[,2]-time[,1])[c[,1]==1]
event.aux<- c[,3][c[,1]==1]

time23 <- inla.surv(time=c(rep(NA,n),rep(NA,n),time.aux), event=c(rep(NA,n),rep(NA,n),event.aux))

#Three intercepts, one for each transition
beta0 <- as.factor(c(rep(1,n), rep(2,n), rep(3,refrac)))

#Covariates

#Age
edad<- edad_i-mean(edad_i)
#Health Area
dpto<- dpto_asignacioncoddesc
dpto<- substr(dpto,1,2)
dpto[substr(dpto,2,2)==":"]<- substr(dpto[substr(dpto,2,2)==":"],1,1)
dpto<- as.numeric(dpto)
#Sex
sexo.aux<-sexo
sexo.aux[sexo.aux=="H"]<- 0
sexo.aux[sexo.aux=="M"]<- 1

#Covariates for transitions F->R and F->D
sex1 <- c(as.numeric(sexo.aux),rep(NA,n),rep(NA,refrac))
sex2 <- c(rep(NA,n),as.numeric(sexo.aux),rep(NA,refrac))
age1 <- c(edad, rep(NA,n),rep(NA,refrac))
age2 <- c(rep(NA,n),edad,rep(NA,refrac))

#Covariates for transition R->D
sex23 <- c(rep(NA,n),rep(NA,n),as.numeric(sexo.aux[c[,1]==1]))
age23 <- c(rep(NA,n),rep(NA,n),edad[c[,1]==1])

dep.aux<- c(dpto,dpto+24)
dep.aux<- c(dep.aux,dpto[c[,1]==1]+48)

detach(data_cohort)

joint.data<- list(sex1=sex1, sex2=sex2, sex23=sex23, age1=age1, age2=age2, age23=age23, beta0=beta0, dep.aux=dep.aux)
joint.data$Y<- list(timeC1,timeC2,time23)

return(joint.data)
}

