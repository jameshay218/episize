#' Mutliple age and titer class final size equation system
#' 
#' Given an attack rate matrix, R0, contact matrix, population proportions and immunity, gives the difference of the final size equation(see Kucharski PLoS Pathogens 2014; Andreasen 2011 and Wallinga 2006). This should return zero if the attack rate matrix is correct.
#' @param A an NxM matrix of attack rates, where N is the number of age classes and M the number of immunity classes
#' @param R0 the disease specific R0 ie. beta/gamma. Note that another parameter will mediate the contact rate
#' @param G the normalised contact rate matrix scaled by population sizes. See \code{\link{setup_C}}
#' @param P NxM matrix of population proportions; number of each age and titre class as proportion of entire population
#' @param Z NxM matrix of immunity. ie. each element is the proportion of infection prevented due to immunity for that age/titre class
#' @return difference between the LHS and RHS of the final size equation
#' @export
simeq <- function(A,R0,G,P,Z){
    f1 <- A - (1-exp(-R0*Z*((G*P)%*%A)))
    f1
}

#' Normalised contact rate matrix
#' 
#' Given a matrix of contact frequency between age classes and a matrix of true population sizes, gives the scaled, normalised contact rate matrix for each age/titre class.
#' @param C1 the non-normalised contact matrix of contact frequencies between each age class
#' @param Ns the matrix of population sizes for each age/titre combination (non-normalised)
#' @return a matrix of the normalised contact rates between all age/titre class combinations.
#' @export
#' #' @examples
#' C <- matrix(c(2,.5,.4,0.3),ncol=2,nrow=2)
#' N11 <- 1000
#' N12 <- 1000
#' N21 <- 2000
#' N22 <- 750
#' N <- matrix(c(N11,N21,N12,N22),ncol=2,nrow=2)
#' C1 <- setup_C(C,N)
setup_C <- function(C1, Ns){
    Ntiter <- ncol(Ns)
    Nage <- nrow(Ns)
    
    M <-  kron(C1,ones(Ntiter,Ntiter)) #' Non-normalised contact matrix scaled for each titre class

    propns <- Ns/rowSums(Ns) #' Age/titre population size as proportion of age population size
    propns_1 <- repmat(propns,Ntiter,Nage) #' Expand to give element for each age/titre combination

    propns_2 <- kron(repmat(t(t(rowSums(Ns)/sum(Ns))),1,Ntiter),ones(Ntiter,Nage)) #' Age population size as proportion of total population size

    C <- M*propns_1/propns_2 #' Generate scaled contact rates for age and titer groups

    return(C)
}

#' Epidemic Final Size Calculation
#' 
#' Calculates the final size of an epidemic given 2-dimensional population categorisation eg. age and immunity class. Note that this uses the final size calculation similar to that in Kucharski et al. 2014 PLoS Pathogens.
#' @param C1 the non-normalised contact matrix of contact frequencies between each age class
#' @param R0 the disease specific R0 ie. beta/gamma. Note that another parameter will mediate the contact rate
#' @param Ns the matrix of population sizes for each age/titre combination (non-normalised) (ie. rows = ages, cols = immunity classes)
#' @param alphas a vector of values between 0 and 1 matching the number of immunity classes
#' @return an NxM matrix of attack rates (ie. proportion of susceptibles becoming infected)
#' @seealso \code{\link{epi_ode_size}}
#' @export
epi_final_size <- function(C1, R0, Ns, alphas){
    Ntiter <- ncol(Ns)
    Nage <- nrow(Ns)

    C <- setup_C(C1, Ns)
    
    propns_3 <- as.numeric(t(Ns/sum(Ns))) #' Generate a matrix of Pai
    propns_3 <- repmat(propns_3,Nage*Ntiter,1)

    A0 <- rand(Ntiter*Nage,1) #' Starting seeds for nleqslv

    rep_alphas <- t(repmat(alphas,1,Nage)) #' Format alphas like A0, as each group will have its incidence reduced

    #' Run optimiser to find attack rates for age/titer groups
    final <- array(nleqslv(A0,simeq,G=C,Z=rep_alphas,P=propns_3,R0=R0,control=list(xtol=1e-15,ftol=1e-15,btol=1e-15,maxit=1000))$x)

    #' Format as matrix
    final <- matrix(final,ncol=Ntiter,nrow=Nage,byrow=T)

    return(final)
}

#' General SIR ode
#' 
#' Generic SIR ode function (for use in deSolve), taking an arbritrary number of populations, ages and immunity classes
#' @param t current time, as for ode function
#' @param vector of compartment sizes eg. S11, I11, R11, S12, I12, R12 etc...
#' @param pars vector of R0, 1/gamma and alphas (vector of immunity conferred for each titre class)
#' @param C the normalised contact matrix of contact frequencies between each age and titre classes
#' @param Nage number of age classes
#' @param Ntiter number of titre classes
#' @return a list of compartment size changes, as required by deSolve
#' @export
general_sir <- function(t,y, pars, C, Nage,Ntiter){
  R0 <- pars[1]
  Tg <- pars[length(pars)]
  alphas <- pars[2:(length(pars)-1)]
    
  sir <- matrix(y,ncol=Nage*Ntiter,nrow=3)

  dS <- -((R0/Tg)*alphas*sir[1,]*(sir[2,]%*%t(C))/sum(sir))
  dR <- sir[2,]/Tg
  dI <- -dS - dR
    
  tmp <- as.vector(rbind(dS,dI,dR))
  
  return(list(c(tmp)))
}

#' Epidemic Final Size Calculation ODE
#' 
#' Calculates the final size of an epidemic given 2-dimensional population categorisation eg. age and immunity class using an SIR model
#' @param C1 the non-normalised contact matrix of contact frequencies between each age class
#' @param R0 the disease specific R0 (transmission rate). Note that another parameter will mediate the contact rate
#' @param Tg number of days spent infectious (ie. 1/gamma)
#' @param Ns the matrix of population sizes for each age/titre combination (non-normalised) (ie. rows = ages, cols = immunity classes)
#' @param alphas a vector of values between 0 and 1 matching the number of immunity classes
#' @return an NxM matrix of attack rates (ie. proportion of susceptibles becoming infected)
#' @seealso \code{\link{epi_final_size}}
#' @export
epi_ode_size <- function(C1, R0, Tg, Ns, alphas){
    C <- setup_C(C1, Ns)

    Nage <- nrow(Ns)
    Ntiter <- ncol(Ns)
    
    long_Ns <- as.numeric(t(Ns))
    start <- NULL
    start[1] <- long_Ns[1] - 1
    start[2] <- 1
    start[3] <- 0
    index <- 4
    for(i in 2:length(long_Ns)){
        start[index] <- long_Ns[i]
        index <- index + 1
        start[index] <- 0
        index <- index + 1
        start[index] <- 0
        index <- index + 1
    }
    y <- ode(y=start,t=seq(1,365,by=1),func=general_sir, parms=c(R0,alphas,Tg),C=C,Nage=Nage,Ntiter=Ntiter)

    A <- NULL
    for(i in 1:((ncol(y)-1)/3)){
        A[i] <- y[nrow(y),(i-1)*3+4]/y[1,(i-1)*3+2]
    }
    A <- matrix(A,nrow=Nage,ncol=Ntiter,byrow=T)
    
    return(A)
}

#' 3 age class SIR ode
#' 
#' 3 age class SIR ode function (for use in deSolve)
#' @param t current time, as for ode function
#' @param vector of compartment sizes eg. S1, I1, R1, S2, I2, R2 etc...
#' @param pars vector of R0, alpha (immunity) and duration of infectious period
#' @param C the normalised contact matrix of contact frequencies between each age classes
#' @return a list of compartment size changes, as required by deSolve
#' @seealso \code{\link{sir}}, \code{\link{sir_2}}
#' @export
sir_3 <- function(t,y, pars, C){
  S1 <- y[1]
  I1 <- y[2]
  R1 <- y[3]
  
  S2 <- y[4]
  I2 <- y[5]
  R2 <- y[6]
  
  S3 <- y[7]
  I3 <- y[8]
  R3 <- y[9]
  
  N1 <- S1 + I1 + R1
  N2 <- S2 + I2 + R2
  N3 <- S3 + I3 + R3
  N <- N1+N2+N3
  print(N)
  R0 <- pars[1]
  alpha <- pars[2]
  durI <- pars[3]
  
  dS1 <- -(R0/durI)*alpha*S1*(C[1,1]*I1/(N) + C[1,2]*I2/(N) + C[1,3]*I3/N)
  dI1 <- (R0/durI)*alpha*S1*(C[1,1]*I1/(N) + C[1,2]*I2/(N) + C[1,3]*I3/N) - I1/durI
  dR1 <- I1/durI
  
  dS2 <- -(R0/durI)*alpha*S2*(C[2,1]*I1/(N) + C[2,2]*I2/(N) + C[2,3]*I3/N)
  dI2 <- (R0/durI)*alpha*S2*(C[2,1]*I1/(N) + C[2,2]*I2/(N) + C[2,3]*I3/N) - I2/durI
  dR2 <- I2/durI
  
  dS3 <- -(R0/durI)*alpha*S3*(C[3,1]*I1/(N) + C[3,2]*I2/(N) + C[3,3]*I3/N)
  dI3 <- (R0/durI)*alpha*S3*(C[3,1]*I1/(N) + C[3,2]*I2/(N) + C[3,3]*I3/N) - I3/durI
  dR3 <- I3/durI
  
  
  return(list(c(dS1,dI1,dR1,dS2,dI2,dR2,dS3,dI3,dR3)))
}

#' Single population SIR case - ODE model
#' 
#' 1 age class SIR ode function (for use in deSolve)
#' @param t current time, as for ode function
#' @param vector of compartment sizes ie. S, I and R
#' @param pars vector of R0 and duration of infectious period
#' @param C the contact rate
#' @return a list of compartment size changes, as required by deSolve
#' @seealso \code{\link{sir_3}}, \code{\link{sir_2}}
#' @export
sir <- function(t,y,pars,C){
  S <- y[1]
  I <- y[2]
  R <- y[3]
  N <- S + I + R
  
  R0 <- pars[1]
  durI <- pars[2]
  
  dS <- -(R0/durI)*S*I*C/N
  dI <- (R0/durI)*S*I*C/N - I/durI
  dR <- I/durI
  return(list(c(dS,dI,dR)))
  
}


#' Two population SIR case - ODE model
#' 
#' 2 age class SIR ode function (for use in deSolve)
#' @param t current time, as for ode function
#' @param vector of compartment sizes ie. S1, I1, R1, S2 etc...
#' @param pars vector of R0, alpha (immunity) and duration of infectious period
#' @param C the 2x2 contact rate matrix
#' @return a list of compartment size changes, as required by deSolve
#' @seealso \code{\link{sir_3}}, \code{\link{sir}}
#' @export
sir_2 <- function(t,y, pars, C){
  S1 <- y[1]
  I1 <- y[2]
  R1 <- y[3]
  
  S2 <- y[4]
  I2 <- y[5]
  R2 <- y[6]
  
  N1 <- S1 + I1 + R1
  N2 <- S2 + I2 + R2
  
  R0 <- pars[1]
  alpha <- pars[2]
  durI <- pars[3]
  
  dS1 <- -(R0/durI)*alpha*S1*(C[1,1]*I1/(N1+N2) + C[1,2]*I2/(N1+N2))
  dI1 <- (R0/durI)*alpha*S1*(C[1,1]*I1/(N1+N2) + C[1,2]*I2/(N1+N2)) - I1/durI
  dR1 <- I1/durI
  
  dS2 <- -(R0/durI)*alpha*S2*(C[2,1]*I1/(N1+N2) + C[2,2]*I2/(N1+N2))
  dI2 <- (R0/durI)*alpha*S2*(C[2,1]*I1/(N1+N2) + C[2,2]*I2/(N1+N2)) - I2/durI
  dR2 <- I2/durI
  
  return(list(c(dS1,dI1,dR1,dS2,dI2,dR2)))
}

#' Two age/two titre class SIR ODE function
#' 
#' Calculates the final size of an epidemic given 2-dimensional population categorisation eg. age and immunity class using an SIR model. Only takes 2 ages and 2 titres.
#' @param t current time, as for ode function
#' @param vector of compartment sizes ie. S11, I11, R11, S12 etc...
#' @param pars vector of R0, alpha1, alpha2 (immunity of titre class 1 and 2) and duration of infectious period
#' @param C the contact rate matrix for contact rates between each age/titre class (4x4 matrix)
#' @return a list of compartment size changes, as required by deSolve
#' @seealso \link{\code{epi_ode_size}}
#' @export
sir_22 <- function(t,y, pars, C){
  S11 <- y[1]
  I11 <- y[2]
  R11 <- y[3]
  
  S12 <- y[4]
  I12 <- y[5]
  R12 <- y[6]
  
  S21 <- y[7]
  I21 <- y[8]
  R21 <- y[9]
  
  S22 <- y[10]
  I22 <- y[11]
  R22 <- y[12]
  
  N11 <- S11 + I11 + R11
  N21 <- S21 + I21 + R21
  N12 <- S12 + I12 + R12
  N22 <- S22 + I22 + R22
  N <- N11 + N21 + N12 + N22
  N1 <- N11 + N12
  N2 <- N21 + N22
  R0 <- pars[1]
  alpha1 <- pars[2]
  alpha2 <- pars[3]
  durI <- pars[4]

  dS11 <- -(R0/durI)*alpha1*S11*((C[1,1]*I11 + C[1,2]*I12 + C[1,3]*I21 + C[1,4]*I22)/N)
  dI11 <- (R0/durI)*alpha1*S11*((C[1,1]*I11 + C[1,2]*I12 + C[1,3]*I21 + C[1,4]*I22)/N) - I11/durI
  dR11 <- I11/durI
  
  dS12 <- -(R0/durI)*alpha2*S12*((C[2,1]*I11 + C[2,2]*I12 + C[2,3]*I21 + C[2,4]*I22)/N)
  dI12 <- (R0/durI)*alpha2*S12*((C[2,1]*I11 + C[2,2]*I12 + C[2,3]*I21 + C[2,4]*I22)/N)- I12/durI
  dR12 <- I12/durI
  
  
  dS21 <- -(R0/durI)*alpha1*S21*((C[3,1]*I11 + C[3,2]*I12 + C[3,3]*I21 + C[3,4]*I22)/N)
  dI21 <- (R0/durI)*alpha1*S21*((C[3,1]*I11 + C[3,2]*I12 + C[3,3]*I21 + C[3,4]*I22)/N) - I21/durI
  dR21 <- I21/durI
  
  dS22 <- -(R0/durI)*alpha2*S22*((C[4,1]*I11 + C[4,2]*I12 + C[4,3]*I21 + C[4,4]*I22)/N)
  dI22 <- (R0/durI)*alpha2*S22*((C[4,1]*I11 + C[4,2]*I12 + C[4,3]*I21 + C[4,4]*I22)/N) - I22/durI
  dR22 <- I22/durI

  return(list(c(dS11,dI11,dR11,dS12,dI12,dR12, dS21,dI21,dR21,dS22,dI22,dR22)))
}
