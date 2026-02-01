adsasi_0d = function(simfun,tar_power=0.9,...,nsims=5000, verbose=FALSE, impNN=Inf, capNN=2000, initiation = TRUE, savegraphs = FALSE, keepsims = FALSE) 
 {
  # Initializing some variables
  par_bak <- par()[c("mfrow","cex","cex.lab","cex.axis")] # graphical parameter backup
  tar_NN = round(exp(seq(log(10),log(300),length.out=48)))    # tar_NN will be a vector of target sample sizes that will be simulated iteratively
  latest_rootsize = sqrt(30)                                 # this is the square root of the sample size
  latest_logslope=-1                                            # this is the slope nuisance parameter (how fast power drops if sample size is not adequate)
  trials = matrix(c(1,sqrt(1000000),0,1),nrow=2)             # this is where we store the TRUE/FALSE results of our simulations for each tested sample size, 
                                                             #    we are initiating assuming 1 patient is not enough and 1M patients is enough
  colnames(trials) = c("srsampsize","signif")            # srsampsize for "Square Root of SAMPle SIZE", signif for whether a discovery is made (even 
                                                             #    if it is not strictly statistical significance for rejection of a well-defined null)
  se = NA                                                    # standard error of the square root of the target sample size using the Hessian matrix
  batch = 0                                                  # just an iteration of simulation batches to keep count
  
  intercept_target = qnorm(tar_power)                        # given the equation we use (see paper), we have an intercept term where we have tar_power = pnorm(intercept_target)

  # Because we use a slightly modified probit regression we have to write the likelihood by hand
  # xx is a scalar vector with two components logslope and rootsize. For any pair of values
  #    in xx, the partial likelihood for "successful" and "failed" trials is computed with 
  #    the normal CDF. Then we sum the partial log-likelihoods to get the total log-likelihood. 
  loglik = function(xx) {sum(
                              pnorm(intercept_target+exp(xx["logslope"])*(trials[!!trials[,"signif"],"srsampsize"]-xx["rootsize"]),lower.tail=TRUE,log.p=TRUE)
                             ,pnorm(intercept_target+exp(xx["logslope"])*(trials[ !trials[,"signif"],"srsampsize"]-xx["rootsize"]),lower.tail=FALSE,log.p=TRUE)
                             )}
  dloglik = function(xx) # gradient
   {
    # see paper for the expression of the gradient for slope and size coefficients
    ss_slope = grep("slope",names(xx))
    ss_size = grep("size",names(xx))
    ordered_vector = xx[c(ss_slope,ss_size)] 
    inner = intercept_target+exp(ordered_vector[1])*(trials[,"srsampsize"]-ordered_vector[2]) # 
    signif_sign = c(-1,+1)[1+trials[,"signif"]] # whether upper (+1) or lower (-1) tail of cdf is taken, depending on significant trial or not ; therefore, also sign of dnorm for derivative
    pnorm_inner = sapply(1:nrow(trials),function(yy){pnorm(inner[yy],lower.tail=as.logical(trials[yy,"signif"]))})
    dnorm_inner = dnorm(inner)
    gradient_slope = (inner-intercept_target)*signif_sign*dnorm_inner/pnorm_inner # gradient wrt to slope
    gradient_size = -exp(ordered_vector[1])*signif_sign*dnorm_inner/pnorm_inner # gradient wrt to size 
    gradient = c(sum(gradient_slope),sum(gradient_size))
    names(gradient) = names(ordered_vector)
    gradient = gradient[names(xx)]
    gradient
    }
    # a bit of code to test the gradient # trials = cbind(srsampsize=c(3,8,10,15,19),signif=c(0,0,1,1,1)) ; xx1 = c(logslope=-3,rootsize=10) ; output = matrix(0,nrow=2,ncol=2) ; for(ii in 1:length(xx1)) { xx2=xx1 ; xx2[ii]=xx1[ii]+.00000001 ; output[ii,1]=((loglik(xx2)-loglik(xx1))/(xx2[ii]-xx1[ii]))[1] ; output[ii,2]=dloglik(xx1)[ii]} ; print(output) ; print(output[,1]/output[,2])
  ddloglik = function(xx) # xx is the same vector, now we compute the Hessian (matrix with second derivatives and derivatives of gradient wrt other parameters)
   {
    # see paper for the expression of the Hessian for slope and size coefficients
    ss_slope = grep("slope",names(xx))
    ss_size = grep("size",names(xx))
    ordered_vector = xx[c(ss_slope,ss_size)] 
    inner = intercept_target+exp(ordered_vector[1])*(trials[,"srsampsize"]-ordered_vector[2]) # 
    signif_sign = c(-1,+1)[1+trials[,"signif"]] # whether upper (+1) or lower (-1) tail of cdf is taken, depending on significant trial or not ; therefore, also sign of dnorm for derivative
    pnorm_inner = sapply(1:nrow(trials),function(yy){pnorm(inner[yy],lower.tail=as.logical(trials[yy,"signif"]))})
    dnorm_inner = dnorm(inner)
    bigfraction = (-inner*dnorm_inner*pnorm_inner-signif_sign*dnorm_inner^2)/(pnorm_inner^2)
    
    hessian_size_size = signif_sign*exp(2*ordered_vector[1])*bigfraction # 
    hessian_size_slope = -signif_sign*exp(ordered_vector[1])*(dnorm_inner/pnorm_inner+(inner-intercept_target)*bigfraction) # 
    hessian_slope_slope = signif_sign*(inner-intercept_target)*(dnorm_inner/pnorm_inner+(inner-intercept_target)*bigfraction) # common term for the region of the Hessian with two slope coeffs
    hessian = matrix(c(sum(hessian_slope_slope),sum(hessian_size_slope),sum(hessian_size_slope),sum(hessian_size_size)),nrow=2,ncol=2)
    dimnames(hessian) = list(rows=names(ordered_vector),cols=names(ordered_vector))
    hessian = hessian[names(xx),names(xx)]
    hessian
    }
    # a bit of code to test the Hessian # output = cbind(xx1,xx2) ; for(ii in names(xx1)) { yy = ii ; xx2=xx1 ; xx2[yy]=xx1[yy]+.00000001 ; output[yy,1]=((dloglik(xx2)-dloglik(xx1))/(xx2[yy]-xx1[yy]))[2] ; output[yy,2]=ddloglik(xx1)[2,yy]} ; print(output) ; print(output[,1]/output[,2]) # repeat with instead of 7 to check another row
    
  # Here is the main loop
  while(  nrow(trials)<nsims                                         # exit if enough simulations made
        & !(nrow(trials)>500 & !is.na(latest_rootsize) & latest_rootsize>sqrt(impNN))  # exit if after 500+ simulated trials the best size seems to be >impNN
        & !(!is.na(latest_rootsize) & latest_rootsize>314 & mean(trials[,"srsampsize"]==sqrt(capNN))>10))  # exit in any case if answer seems >100k patients and the simulator is hitting capNN (likely a user error leading to lots of wasted compute)
   {
    if(nrow(trials)>200&&is.null(dim(initiation))&&initiation==FALSE) { trials=trials[-(1:150),] ; initiation = TRUE } # kick out first iterations, deactivate the logical switch
    
    # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # 
    # Running next simulations
    # # # # # # # # # # # # # # # # # # 
    args = list(NN=0,...)  # preparing arguments for vectorized simfun call
    args = lapply(1:length(tar_NN),function(xx){yy=args;yy[[1]]=tar_NN[xx];yy}) # this is a list of lists, each list has the arguments to be passed to simfun
    
    simulations = unlist(lapply(args,function(xx){do.call(simfun,xx)})) # running the simulations, normally getting a vector of logicals. Note that this calls what you have defined to be simfun when calling adsasi and passes whatever arguments you have specified, except for sample sizes (values in tar_NN are used)
    trials = rbind(trials,cbind(srsampsize=sqrt(tar_NN),simulations))  # collating the latest simulations with the ones before
    
    # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # 
    # Fitting model to simulation results
    # # # # # # # # # # # # # # # # # # # # 
    initial_vector = c(logslope=-5,latest_rootsize) # if logslope is too small or too big the gradient is invisible
    names(initial_vector) = c("logslope","rootsize")
    initial_value = loglik(initial_vector)    # log-likelihood for the initial values passed to the minimizer 
    if(initial_value!=Inf&initial_value!=-Inf&!is.na(initial_value))    # if the initial value is defined, the minimizer can get to work
     {
      try(mle <- optim(initial_vector, fn=loglik, gr=dloglik, method="BFGS", hessian = FALSE, control=list(fnscale=-1)), silent = TRUE) # using the default optimizer in R, a bit slow but the main compute use is for within-simulation inference anyway
      # control=list(fnscale=-1)) to maximize instead of minimizing
      # here we use try() in case it fails, gets stuck etc. in which case we keep previous values
      latest_logslope = mle$par["logslope"]
      latest_rootsize = mle$par["rootsize"]

      covmat = diag(length(mle$par))
      try(covmat <- solve(-ddloglik(c(latest_logslope,latest_rootsize))), silent=TRUE)             # using the exact Hessian to get a Rao Cramer covariance matrix, will fail from time to time, mostly when not invertible
      
      # Computing confidence intervals using variances for logslope and rootsize
      slope_estimate = exp(latest_logslope)
      slope_confint_lower = exp(latest_logslope-1.96*sqrt(covmat[1,1]))
      slope_confint_higher = exp(latest_logslope+1.96*sqrt(covmat[1,1]))
      size_estimate = unname(latest_rootsize^2)
      se = sqrt(covmat[2,2]) ; if(se==0) { se=.2*latest_rootsize }
      size_confint_lower = (latest_rootsize-1.96*se)^2
      size_confint_higher = (latest_rootsize+1.96*se)^2
      
      batch_size = max(50,round(nrow(trials)/10))              # the number of simulations for the new batch is 10% of what's already been done, but at least 50
      if((nrow(trials)+batch_size)>nsims) batch_size = max(2,nsims-nrow(trials)) # max(2,...) is to avoid bugs so one may rarely get nsims+1 runs at the end, depending on nsims
      tar_NN = round(seq(latest_rootsize-se,latest_rootsize+se,length.out=batch_size)^2)
      tar_NN = pmin(pmax(2,tar_NN),round(runif(length(tar_NN),capNN*0.5,capNN)))              # drawing new tar_NNs to try, applying limits : minimum N=4, maximum=capNN
      } else {                                                 # if optimizer not launchable, keep the same estimates + some noise, keep target trials (N & optival)
              latest_logslope = latest_logslope+rnorm(1)
              latest_rootsize = latest_rootsize+rnorm(1)
              warning(paste(sep="","Adsasi wrapper log-likelihood fit fail after ",nrow(trials)," simulations."))  # returning a warning, this is usually bad if more than once/twice per run
              }       
        
    if(all(c(trials[,1]^2,tar_NN)==tar_NN[1])) tar_NN = tar_NN*rep(c(.5,1),c(round(length(tar_NN)),length(tar_NN)-round(length(tar_NN))))  # we divide half of values by 2 if all trials have the same size
    if(verbose) { cat("\n") ; print(tar_NN); cat("\n") ; print(se) ; cat("\n") }    # another verbose descriptions of how the run is going
    
    # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # 
    # Plotting
    # # # # # # # # # # # # # # # # # # 
    if(savegraphs!=FALSE) { batch=batch+1 ; pngname = paste(sep="","adsasi",gsub(":","-",as.character(Sys.time()),fixed=TRUE)," iteration ",paste(rep(0,3-nchar(batch)),collapse=""),batch,".png") ; if(is.character(savegraphs)){pngname=gsub("adsasi",savegraphs,pngname)} ; png(pngname,width=400*2,height=400) } # opening a device to save graphs
    #layout(mat = matrix(c(1,1,1,1,2,3),nrow=2), heights = c(1,1), widths = c(1,1,1), respect =TRUE)
    par(mfrow=c(1,2),cex=1.35,cex.lab=1,cex.axis=1)   # light aesthetic setup, with silent backup per standard R usage
    on.exit(par(par_bak)) # restoring on exit
    plot(
          trials[,1]^2                                        # sample sizes as ordinates, the abscissae are by default 1:nrow(trials) and that's what we want (to plot them in order)
         ,xlim=c(0,nsims*1.1)
         ,ylim=quantile(trials[,1]^2,c(.1,.9)) * c(.75,1.33)  # zoomed, not showing everything
         ,main=c("Trace of sample size exploration",paste0("Current estimate ",round(size_estimate,1)," = (",round(latest_rootsize,1),"+-",round(se,2),")^2")) # writing the latest estimate and its error
         ,xlab="Simulation#",ylab="Sample size (exact)",type="p"
         ,col=paste(sep="",c(rgb(1,.5,0),rgb(0,.6,0.8)),"")[1+trials[,2]],pch=3,cex=.5+.1*(1-trials[,2]) # successes in blue, failures in orange ; so higher ordinates are bluer, lower are oranger
         )

    points(nrow(trials)+1:length(tar_NN),tar_NN,col="#55555588",pch=1,cex=.8)   # drawing current batch on the graph, result unknown so in grey
    plot(                                                                     # second graph (right panel)
          trials[,1]^2+rnorm(nrow(trials),0,.25)                              # sample size
         ,trials[,2] + (trials[,2]-.5)*2*-rexp(nrow(trials),20)               # success or failure, with some jittering to appreciate density
         ,xlim=latest_rootsize^2 * c(.75,1.33)                                # slightly zoomed, not showing everything
         ,ylim=0:1,main=c("Probit regression",paste0("New estimate ",round(size_estimate,1),", slope ",round(exp(latest_logslope),3))),xlab="Sample size (jittered)",ylab="Expected power"
         ,type="p",col=paste(sep="",c(rgb(1,.5,0),rgb(0,.6,0.8)),"")[1+trials[,2]],pch=3,cex=.5
         )
    lines(latest_rootsize^2*c(.75,1,1),c(tar_power,tar_power,0),lty=2)
    xx = seq(latest_rootsize^2*.75,latest_rootsize^2*1.33,length.out=1000)
    lines(xx,pnorm(qnorm(tar_power)+exp(latest_logslope)*(sqrt(xx)-latest_rootsize)))

    if(savegraphs!=FALSE) dev.off()                                                  # closing device and saving image if graphs are to be saved
    }                                                                         # end of loop
  
  size_estimate = ceiling(size_estimate)
  if(nrow(trials)<nsims) warning("Halted for inefficiency.")
  if( size_estimate > 4*capNN ) { warning("Required sample size seems much larger than user-provided size cap, returning +Inf.") ; size_estimate=Inf }  # there is little accuracy when trying to extrapolate that far 
  output = list(size_estimate=size_estimate)
  if(keepsims) 
   {
    extra_output = list( trials = trials, confint=unname(c(size_confint_lower,size_confint_higher)) )
    output = c(output,extra_output)
    }
  output
  }
