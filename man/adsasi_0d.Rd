\name{adsasi_0d}
\alias{adsasi_0d}
\title{Adaptive Sample Size Finder For Fixed Designs}
\usage{
adsasi_0d(
  simfun,
  tar_power = 0.9,
  ...,
  nsims = 5000,
  verbose = FALSE,
  impNN = Inf,
  capNN = 2000,
  initiation = TRUE,
  savegraphs = FALSE,
  keepsims = FALSE
)
}
\arguments{
\item{simfun}{(function) The user-supplied function that describes the clinical trial scenario (or similar experiment) that needs to be explored. Must have as named arguments a sample size (named \code{NN}) and an arbitrary number of design parameters. Must return a boolean indicating whether the trial is successful or not, after performing any required computations (regressions, bootstraps) as written by the user, and never return \code{NA}.}

\item{tar_power}{(single number between 0 and 1) Target power (or more broadly, probability of success). \code{adsasi_0d} will seek regions where \code{simfun} returns \code{TRUE} with a frequency of \code{tar_power}, assuming that higher sample size equals higher probability of success.}

\item{...}{Additional named arguments to be passed to \code{simfun}. Some of these arguments can be functions themselves (e.g. for trying different analysis models). Any \code{simfun} argument without a default value must be specified here.}

\item{nsims}{(single number) Number of simulations to be run. After initialization, simulations are run in batches of 10\% of the number of existing simulations, until \code{nsims} is reached.}

\item{verbose}{(boolean) Whether to print extra diagnostics messages throughout the run.}

\item{impNN}{(single number, or infinity) Sample size that is considered impossible (either computationnally, or logistically). The simulator will exit if, after 500+ simulations, it looks like the best value is above this. In practice, is mostly useful to avoid expensive computations in situations where \code{simfun} is not written well or is prohibitively long to compute for large sample sizes.}

\item{capNN}{(single number, or infinity) Maximum sample size that will be simulated. Also mostly useful to avoid expensive computations. Values between \code{capNN} and \code{impNN} will be extrapolations of unclear validity, so if it looks like the answer is really above \code{capNN}, try running the wrapper again with a higher \code{capNN}.}

\item{initiation}{(boolean, or numeric matrix) Either a boolean indicating whether or not to keep the first 150 simulations for the relationship inference (those tend to be far from \code{tar_power}), or a matrix with simulation results from a previous run which the user wants enrich with more simulations (formatted exactly as produced by \code{adsasi_0d} with the same \code{simfun}). See \code{keepsims} and Note below for how to store and retrieve this data.}

\item{savegraphs}{(boolean or string) Whether to save graphs on drive (vs. showing them in the console). If string, is interpreted as a typical name to be used (several graphs will be drawn, with iteration number, timestamp and .png file extension appended). The string can contain a filepath, but folders must already exist (e.g. with \code{dir.create()} from \code{base}, if automated).}

\item{keepsims}{(boolean) Whether to keep the simulations sizes and individual outcomes in the output. See Note for format details.}
}
\value{
(list) A list with one (by default) element named \code{size_estimate} with the sample size to obtain probability of success equal to \code{tar_power}. If \code{keepsims=TRUE}, additional elements with fits and simulation results (see Note).  
}
\description{
This function empirically finds the relationship between sample size and power, for a given experiental simulation scenario supplied by the user in the form of a function (most typically a clinical trial, but any experiment whose success rate increases with the number of observations can be processed). \code{adsasi_0d} will try different sample sizes and progressively zoom on the ones where power is nominal. Power is understood in a broad sense here, as a probability of success of the experiment rather than a strict statistical power.
}
\note{
Standard error of estimated sample size is shown on graphs but not saved by default anywhere. 

With \code{keepsims=TRUE}, additional elements will be returned : \code{trials} for all the used simulations, \code{confint} for 95\% confidence interval of \code{size_estimate}. The confidence interval is obtained in square root form and is assymetric ; the reverse calculation needs to be done to extract the standard error of its root (take the roots of the bounds & divide their distance by 2*1.96). 

The \code{trials} element can be fed back into \code{adsasi_0d} using argument \code{initiation=trials} to resume with the same simulations as before. Note that the second \code{adsasi_0d} call will need to use the same \code{simfun} and its arguments for this to make sense.
}
\examples{
# First, the user defines a function for their target situation. In this simple example, a 2-sample
# t-test with unequal allocation. Note the syntax to avoid returning NAs. 
simulate_unequal_t_test = function(NN=20,ratio_n1_NN=0.5,delta=1)
 {
  n1 = round(ratio_n1_NN*NN) ; n2 = NN-n1
  yy1 = rnorm(n1) ; yy2 = rnorm(n2,delta)
  pp=NA ; try(pp <- t.test(yy1,yy2)$p.value,silent=TRUE)
  !is.na(pp) & pp<0.05
  }
simulate_unequal_t_test()
# Now we empirically find the relationship between sample size and the parameter of interest. 
# Note that we can change the simfun parameters directly from the adsasi_0d call. 
# nsims should generally be much higher than in this fast-running example (>5000). 
adsasi_0d(simulate_unequal_t_test,delta=1.25,nsims=200)
}
\keyword{optimization}
\keyword{sample}
\keyword{simulation}
\keyword{size}
