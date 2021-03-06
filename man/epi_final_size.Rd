% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/final_size.R
\name{epi_final_size}
\alias{epi_final_size}
\title{Epidemic Final Size Calculation}
\usage{
epi_final_size(C1, R0, Ns, alphas)
}
\arguments{
\item{C1}{the non-normalised contact matrix of contact frequencies between each age class}

\item{R0}{the disease specific R0 ie. beta and gamma. Note that another parameter will mediate the contact rate}

\item{Ns}{the matrix of population sizes for each age or titre combination (non-normalised) (ie. rows = ages, cols = immunity classes)}

\item{alphas}{a vector of values between 0 and 1 matching the number of immunity classes}
}
\value{
an NxM matrix of attack rates (ie. proportion of susceptibles becoming infected)
}
\description{
Calculates the final size of an epidemic given 2-dimensional population categorisation eg. age and immunity class. Note that this uses the final size calculation similar to that in Kucharski et al. 2014 PLoS Pathogens.
}
\seealso{
\code{\link{epi_ode_size}}
}

