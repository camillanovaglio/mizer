% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setReproduction.R
\name{setReproduction}
\alias{setReproduction}
\alias{getMaturityProportion}
\alias{getReproductionProportion}
\title{Set reproduction parameters}
\usage{
setReproduction(params, maturity = NULL, repro_prop = NULL, RDD = NULL, ...)

getMaturityProportion(params)

getReproductionProportion(params)
}
\arguments{
\item{params}{A MizerParams object}

\item{maturity}{Optional. An array (species x size) that holds the proportion
of individuals of each species at size that are mature. If not supplied, a
default is set as described in the section "Setting reproduction".}

\item{repro_prop}{Optional. An array (species x size) that holds the
proportion of consumed energy that a mature individual allocates to
reproduction for each species at size. If not supplied, a default is set as
described in the section "Setting reproduction".}

\item{RDD}{The name of the function calculating the density-dependent
reproduction rate from the density-independent rate. Defaults to
"\code{\link[=BevertonHoltRDD]{BevertonHoltRDD()}}".}

\item{...}{Unused}
}
\value{
The updated MizerParams object. Because of the way the R language
works, \code{setReproduction()} does not make the changes to the params object
that you pass to it but instead returns a new params object. So to affect
the change you call the function in the form
\code{params <- setReproduction(params, ...)}.
}
\description{
Sets the proportion of the total energy available for reproduction and growth
that is invested into reproduction as a function of the size of the
individual and sets the reproductive efficiency.
}
\section{Setting reproduction}{


\subsection{Investment}{
For each species and at each size, the proportion of the available energy
that is invested into reproduction is the product of two factors: the
proportion \code{maturity} of individuals that are mature and the proportion
\code{repro_prop} of the energy available to a mature individual that is
invested into reproduction.

If the \code{maturity} argument is not supplied, then it is set to a sigmoidal
maturity ogive that changes from 0 to 1 at around the maturity size:
\deqn{{\tt maturity}(w) = \left[1+\left(\frac{w}{w_{mat}}\right)^{-U}\right]^{-1}.}{
  maturity(w) = [1+(w/w_mat)^(-U)]^(-1)}
(To avoid clutter, we are not showing the species index in the equations.)
The maturity weights are taken from the \code{w_mat} column of the
species_params data frame. Any missing maturity weights are set to 1/4 of the
asymptotic weight in the \code{w_inf} column.
The exponent \eqn{U} determines the steepness of the maturity ogive. By
default it is chosen as \eqn{U = 10}, however this can be overridden by
including a column \code{w_mat25} in the species parameter dataframe that
specifies the weight at which 25\% of individuals are mature, which sets
\eqn{U = \log(3) / \log(w_{mat} / w_{25}).}{U = log(3) / log(w_mat / w_25).}

The sigmoidal function given above would strictly reach 1 only asymptotically.
Mizer instead sets the function equal to 1 already at the species'
maximum size, taken from the compulsory \code{w_inf} column in the
\code{species_params} data frame.

If the \code{repro_prop} argument is not supplied, it is set to the
allometric form
\deqn{{\tt repro\_prop}(w) = \left(\frac{w}{w_{inf}}\right)^{m-n}.}{
  repro_prop = (w/w_inf)^(m - n).}
Here \eqn{n} is the scaling exponent of the energy income rate. Hence
the exponent \eqn{m} determines the scaling of the investment into
reproduction for mature individuals. By default it is chosen to be
\eqn{m = 1} so that the rate at which energy is invested into reproduction
scales linearly with the size. This default can be overridden by including a
column \code{m} in the species parameter dataframe. The asymptotic sizes
are taken from the compulsory \code{w_inf} column in the species_params
data frame.
}

\subsection{Efficiency}{
The reproductive efficiency, i.e., the proportion of energy allocated to
reproduction that results in egg biomass, is set from the \code{erepro}
column in the species_params data frame. If that is not provided, the default
is set to 1 (which you will want to override). The offspring biomass divided
by the egg biomass gives the rate of egg production, returned by
\code{\link[=getRDI]{getRDI()}}.
}

\subsection{Density dependence}{
The stock-recruitment relationship is an emergent phenomenon in mizer, with
several sources of density dependence. Firstly, the amount of energy invested
into reproduction depends on the energy income of the spawners, which is
density-dependent due to competition for prey. Secondly, the proportion of
larvae that grow up to recruitment size depends on the larval mortality,
which depends on the density of predators, and on larval growth rate, which
depends on density of prey.

Finally, to encode all the density dependence in the stock-recruitment
relationship that is not already included in the other two sources of density
dependence, mizer puts the the density-independent rate of egg production
through a density-dependence function. The result is returned by
\code{\link[=getRDD]{getRDD()}}. The name of the density-dependence function is
specified by the \code{RDD} argument. The default is the Beverton-Holt
function \code{\link[=BevertonHoltRDD]{BevertonHoltRDD()}}, which requires an \code{R_max} column
in the species_params data frame giving the maximum egg production rate. If
this column does not exist, it is initialised to \code{Inf}, leading to no
density-dependence. Other functions provided by mizer are
\code{\link[=RickerRDD]{RickerRDD()}} and \code{\link[=SheperdRDD]{SheperdRDD()}} and you can easily use
these as models for writing your own functions.
}
}

\examples{
\dontrun{
params <- NS_params
# Change maturity size for species 3
params@species_params$w_mat[3] <- 24
params <- setReproduction(params)
}
}
\seealso{
Other functions for setting parameters: 
\code{\link{setExtMort}()},
\code{\link{setFishing}()},
\code{\link{setInitialValues}()},
\code{\link{setInteraction}()},
\code{\link{setMaxIntakeRate}()},
\code{\link{setMetabolicRate}()},
\code{\link{setPredKernel}()},
\code{\link{setResource}()},
\code{\link{setSearchVolume}()},
\code{\link{species_params}()}
}
\concept{functions for setting parameters}
