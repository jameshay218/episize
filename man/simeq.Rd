% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/final_size.R
\name{simeq}
\alias{simeq}
\title{Mutliple age and titer class final size equation system}
\usage{
simeq(A, R0, G, P, Z)
}
\arguments{
\item{A}{an NxM matrix of attack rates, where N is the number of age classes and M the number of immunity classes}

\item{R0}{the disease specific R0 ie. beta/gamma. Note that another parameter will mediate the contact rate}

\item{G}{the normalised contact rate matrix scaled by population sizes. See \code{\link{setup_C}}}

\item{P}{NxM matrix of population proportions; number of each age and titre class as proportion of entire population}

\item{Z}{NxM matrix of immunity. ie. each element is the proportion of infection prevented due to immunity for that age/titre class}
}
\value{
difference between the LHS and RHS of the final size equation
}
\description{
Given an attack rate matrix, R0, contact matrix, population proportions and immunity, gives the difference of the final size equation(see Kucharski PLoS Pathogens 2014; Andreasen 2011 and Wallinga 2006). This should return zero if the attack rate matrix is correct.
}

