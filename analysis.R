library(deSolve)
my_ode=function(t,state,parms){
  with(as.list(state),{
    dndt=rep(0,length(state))
    dndt[1]=-a*G
    dndt[2]=a*G-K*B
    return(list(dndt))
  })
}
a=0.9
t=seq(0,80,0.1)
G=1000
K=log(2)/5
B=0
init=c(G=G,B=B)
out= ode(y = init, times = t, func = my_ode, parms = NULL)
write.csv(as.data.frame(out), "results/drug_absorption.csv", row.names=FALSE)
head(out) 
plot(out,type="l",xlab="Time (h)",ylab="Drug amount (mg)")

out_df=as.data.frame(out)

max_index=which.max(out_df$B)
time_at_max_concentration <- out_df$time[max_index]

print(time_at_max_concentration)

# Ligand-receptor binding

library(deSolve)
my_ode=function(t,state,parms){
  with(as.list(state),{
    dndt=rep(0,length(state))
    dndt[1]= -kon*R*L +koff*C # receptor
    dndt[2]= kon*R*L -koff*C #complex
    return(list(dndt))
  })
}
R=2000
koff=0.01
kon=0.6
L=2
C=0
t=seq(0,15,0.1)
init=c(R=R,C=C)
out <- ode(y = init, times = t, func = my_ode, parms = NULL)
write.csv(as.data.frame(out), "results/receptor_binding.csv", row.names=FALSE)
plot(out,type="l",xlab="Time (s)",ylab="Number")

#part6
head(out)
in3=out[which(out[,"time"]==3),"C"]



library(deSolve)

############# ODE function #########
my_ode<-function(t,state,parms){
  with(as.list(state),{
    dndt=rep(0,length(state))
    
    dndt[1]= -(beta*S*I)/N+delta*R # dS/dt
    dndt[2]= (beta*S*I)/N -I*gamma # dI/dt
    dndt[3]= I*gamma-delta*R # dR/dt
    return(list(dndt)) # Return 
  })
}
#-------------------------------------
beta=0.4 # Risk of transmitting per contact 
delta=0.005 # Rate of loss of resistance (per day)
N=1000 # Total number of individuals in the population
gamma=0.2
S=995 # Initial number of susceptible
I=5 # Initial number of infected
R=0 # Initial number of resistant
init=c(S=S,I=I,R=R)
t=seq(0,1000,0.1)
out <- ode(y = init, times = t, func = my_ode, parms = NULL)
write.csv(as.data.frame(out), "results/sirs.csv", row.names=FALSE)
Suscept<-out[,2]
Infect<-out[,3]
Resist<-out[,4]
time<-out[,1]
plot(time,Suscept,type="l",xlab="Time (day)",ylab="Number",ylim=c(0,1100),col="pink")
lines(time,Infect,col="yellow")
lines(time,Resist,col="black")
legend("topright",c("Susceptible","Infected","Resistant"),col=c("pink","yellow","black"),lty=1)


library(deSolve)

############# ODE function ####### ##
my_ode<-function(t,state,parms){
  with(as.list(state),{
    dndt=rep(0,length(state))
    
    dndt[1]= -(beta*S*I)/N+delta*R # dS/dt
    dndt[2]=(beta*S*I)/N-omega*E # dE/dt
    dndt[3]=omega*E-gamma*I #dI/dt
    dndt[4]= I*gamma-delta*R # dR/dt
    return(list(dndt)) # Return 
  })
}
#-------------------------------------
beta=0.5 # Risk of transmitting per contact 
delta=0.03 # Rate of loss of resistance (per day)
N=1000# Total number of individuals in the population
gamma=0.1
S=995 # Initial number of susceptible
omega=0.1
I=5 # Initial number of infected
R=0 # Initial number of resistant
E=0 # Initial number of exposed
init=c(S=S,E=E,I=I,R=R)
t=seq(0,200,1)
out <- ode(y = init, times = t, func = my_ode, parms = NULL)
write.csv(as.data.frame(out), "results/seirs.csv", row.names=FALSE)
Suscept<-out[,2]
Infect<-out[,4]
Exposed<-out[,3]
Resist<-out[,5]
time<-out[,1]
plot(time,Suscept,type="l",xlab="Time (day)",ylab="Number",ylim=c(0,1000),col="pink")
lines(time,Infect,col="yellow")
lines(time,Exposed,col="purple")
lines(time,Resist,col="black")
legend("topright",c("Susceptible","Exposed","Infected","Resistant"),col=c("pink","purple","yellow","black"),lty=1)

