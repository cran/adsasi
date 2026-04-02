\name{adsasi_1d}
\alias{adsasi_1d}
\title{Adaptive Sample Size Finder With One Floating Parameter}
\usage{
adsasi_1d(
  simfun,
  tar_power = 0.9,
  ...,
  optivar,
  optiwin = c(min = 0, max = 1),
  optilog = FALSE,
  optiround = FALSE,
  nsims = 5000,
  verbose = FALSE,
  impNN = Inf,
  capNN = 2000,
  initiation = TRUE,
  savegraphs = FALSE,
  keepsims = FALSE,
  n_slope_coefs = 3,
  n_size_coefs = 5
)
}
\arguments{
\item{simfun}{(function) The user-supplied function that describes the clinical trial scenario (or similar experiment) that needs to be explored. Must have as named arguments a sample size (named \code{NN}) and an arbitrary number of design parameters (one of which will be optimized). Must return a boolean indicating whether the trial is successful or not, after performing any required computations (regressions, bootstraps) as written by the user, and never return \code{NA}.}

\item{tar_power}{(single number between 0 and 1) Target power (or more broadly, probability of success). \code{adsasi_1d} will seek regions where \code{simfun} returns \code{TRUE} with a frequency of \code{tar_power}, assuming that higher sample size equals higher probability of success.}

\item{...}{Additional named arguments to be passed to \code{simfun}. Some of these arguments can be functions themselves (e.g. for trying different analysis models). Any \code{simfun} argument without a default value must be specified here.}

\item{optivar}{(single string) Name of the \code{simfun} argument that needs to be optimized.}

\item{optiwin}{(numeric vector of size 2) Bounds of the region to be explored for values of \code{optivar}.}

\item{optilog}{(boolean) Whether \code{optivar} is best explored and drawn in log scale (as in the case of a ratio) or linearly. If for example \code{optiwin=c(0.1,10)} for a ratio, graphs will be drawn with 1 as the middle value if \code{optilog=TRUE} and 5.05 if \code{optilog=FALSE}.}

\item{optiround}{(boolean) Whether \code{optivar} needs to be rounded to the nearest integer to make sense for \code{simfun} (for example, if it is a number of centers in a cluster-randomized trial).}

\item{nsims}{(single number) Number of simulations to be run across all values of \code{optivar}. After initialization, simulations are run in batches of 10\% of existing simulations, until \code{nsims} is reached.}

\item{verbose}{(boolean) Whether to print extra diagnostics messages throughout the run.}

\item{impNN}{(single number, or infinity) Sample size that is considered impossible (either computationnally, or logistically). The simulator will exit if, after 500+ simulations, it looks like the best value is above this. In practice, is mostly useful to avoid expensive computations in situations where \code{simfun} is not written well or is prohibitively long to compute for large sample sizes.}

\item{capNN}{(single number, or infinity) Maximum sample size that will be simulated. Also mostly useful to avoid expensive computations. Initiation will also be done between 10 and \code{capNN/2} patients, so absurdly high values will cause a slow start.}

\item{initiation}{(boolean, or numeric matrix) Either a boolean indicating whether or not to keep the first 150 simulations for the relationship inference (those tend to be far from \code{tar_power}), or a matrix with simulation results from a previous run which the user wants enrich with more simulations (formatted exactly as produced by \code{adsasi_1d} with the same \code{opti*} arguments). See Value and Notes below for how to get this data from the output to be able to reuse it.}

\item{savegraphs}{(boolean or string) Whether to save graphs on drive (vs. showing them in the console). If string, is interpreted as a typical name to be used (several graphs will be drawn, with iteration number, timestamp and .png file extension appended). The string can contain a filepath, but folders must already exist (e.g. with \code{dir.create()} from \code{base}, if automated).}

\item{keepsims}{(boolean or string) Whether to keep simulations and last fit in the returned object, which by default only containe the best value.}

\item{n_slope_coefs}{(single integer) Number of coefficients for the slope polynomial. The slope polynomial tries to model the relationship between \code{optivar} and the loss of power as sample size locally deviates from the (unknown) target.}

\item{n_size_coefs}{(single integer) Number of coefficients for the size polynomial. The size polynomial tries to model the relationship between \code{optivar} and the target sample size. Its shape is the most useful output of the function.}
}
\value{
A \code{list} with 2 numbers in it : minimum sample size, named \code{min_NN}, and corresponding best parameter value, named \code{min_optival}. If \code{keepsims=TRUE}, several other objects will be appended to the list (see Note). 
}
\description{
This function empirically finds the relationship between sample size and a numeric parameter of interest, for a given experiental simulation scenario supplied by the user in the form of a function (most typically a clinical trial, but any experiment whose success rate increases with the number of observations can be processed). \code{adsasi_1d} will search the two-dimensional space empirically (sample size x parameter of interest), favoring exploration of low sample size regions, to find the line where power is nominal. Power is understood in a broad sense here, as a probability of success of the experiment rather than a strict statistical power.
}
\note{
The graph modelling the relationship between parameter value and sample size is generally the most useful output, and is shown but not saved by default.

With \code{keepsims=TRUE}, the function keeps summary simulation results in the returned list, which can, among others, be used to draw the main graph again in a different style (as in Examples). The returned list will have the following extra elements : \code{min_NN} (last sample sizes simulated), \code{min_optival} (corresponding values for the parameter indicated by \code{optivar}, scaled between -1 and 1), \code{trials} (all the used simulations, including a rescaled \code{optival}), \code{abscissae} (natural-scale values for the optimization parameter, for which sample sizes have been computed), \code{slope_natural_estimate_by_optival} (slope variation by \code{optivar} values, see below for plotting), \code{slope_confint_lower_by_optival} (lower bound of confidence interval), \code{slope_confint_higher_by_optival} (higher bound of confidence interval), \code{size_natural_estimate_by_optival} (sample size variation by \code{optivar} values, see below for plotting), \code{size_confint_lower_by_optival} (lower bound of confidence interval), \code{size_confint_higher_by_optival} (higher bound of confidence interval). 

The \code{trials} element can be fed back into \code{adsasi_1d} using argument \code{initiation=x[["trials"]]} (if the previous call was saved in \code{x}) to resume with the same simulations as before. Note that the second \code{adsasi_1d} call will need to use the same \code{simfun}, fixed \code{simfun} arguments, \code{optivar}, \code{optiwin} and \code{optilog} arguments for this to make sense, because values for \code{optivar} stored in \code{trials} are between -1 and +1 and are scaled using these arguments before being passed to \code{simfun} or shown to the user. If one wants to widen the window, the rescaling will need to be done manually.
}
\examples{
# First, the user defines a function for their target situation. In this simple example, a 2-sample
# t-test with unequal allocation. The design parameter of interest will be the ratio of 
# n1 (observations in arm 1) to NN (total sample size). Note the syntax to avoid returning NAs. 
simulate_unequal_t_test = function(NN=20,ratio_n1_NN=0.5,delta=1)
 {
  n1 = round(ratio_n1_NN*NN) ; n2 = NN-n1
  yy1 = rnorm(n1) ; yy2 = rnorm(n2,delta)
  pp=NA ; try(pp <- t.test(yy1,yy2)$p.value,silent=TRUE)
  !is.na(pp) & pp<0.05
  }
simulate_unequal_t_test()
# Now we empirically find the relationship between sample size and the parameter of interest. 
# Note that we can change the delta parameter directly from the adsasi_1d call. 
# nsims should generally be much higher than in this fast-running example (>5000). 
batch=adsasi_1d(simulate_unequal_t_test,delta=1.25,optivar="ratio_n1_NN",nsims=200,keepsims=TRUE)
# Drawing the output in a different style
plot( batch[["abscissae"]],batch[["size_natural_estimate_by_optival"]]
     ,xlab="Optimization parameter",ylab="Estimated sample size",type="o",col="red"
     )
polygon( c(batch[["abscissae"]],rev(batch[["abscissae"]]))
        ,c( batch[["size_confint_higher_by_optival"]]
           ,rev(batch[["size_confint_lower_by_optival"]])
           )
        ,col="#55000088",border=NA
        )
}
\keyword{optimization}
\keyword{sample}
\keyword{simulation}
\keyword{size}
