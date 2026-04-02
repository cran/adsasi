adsasi_1d = function(simfun,tar_power=0.9,...,optivar,optiwin=c(min=0,max=1),optilog=FALSE,optiround=FALSE,nsims=5000, verbose=FALSE, impNN=Inf, capNN=2000, initiation = TRUE, savegraphs = FALSE, keepsims = FALSE, n_slope_coefs=3,n_size_coefs=5) 
 {
  # Initializing some variables
  par_bak <- par()[c("mfrow","cex","cex.lab","cex.axis")]                   # graphical parameter backup
  tar_NN = rep(round(exp(seq(log(10),log(capNN/2),length.out=10))),10)      # tar_NN will be a vector of target sample sizes that will be simulated iteratively
  tar_optival = rep(seq(-1,1,length.out=12)[c(-1,-12)],rep(10,10))          # tar_optival is the vector of optimization parameters to be simulated along with tar_NN. The algorithm searches the -1 to 1 space and converts that into the window provided by the user
  latest_coefs_size = c(10,rep(0,n_size_coefs-1))                           # latest estimates, here this means initial guess of sample size is 10^2=100 for any optival (only first term is nonzero)   
  latest_coefs_slope = c(-10,rep(0,n_slope_coefs-1))                        # latest slopes, here this means slope of exp(-10)=5e-5 for any optival (essentially flat) ; the optimizer can incline a flat slope but has problems flattening a too-steep one  
  trials = cbind(matrix(c(1,sqrt(1000000),0,1),nrow=2)[rep(1:2,20),],seq(-1,1,length.out=42)[c(-1,-42)]) # assuming a 1M sample size yields success whatever the optimization parameter value
  if(!is.null(dim(initiation))) trials = initiation                         # if initiation is a matrix 
  if(any(colnames(trials)=="optival_natural")) trials = trials [,colnames(trials)!="optival_natural"]

  colnames(trials) = c("srsampsize","signif","optival")      # srsampsize for "Square Root of SAMPle SIZE", signif for whether a discovery is made (whether or not strictly a rejection of a well-defined null, the user decides what they want when they write simfun), optival for the optimization parameter value
  se = NA                                                    # standard error of the square root of the target sample size using the Hessian matrix
  batch = 0                                                  # just an iteration of simulation batches to keep count
  
  intercept_target = qnorm(tar_power)                        # given the equation we use (see paper), we have an intercept term where we have tar_power = pnorm(intercept_target)

  # Because we use a slightly modified probit regression we have to write the likelihood, gradient and Hessian by hand
  # xx is a scalar vector with six components for size and slope (polynomial up to ^5). 
  # For any given values of these coefficients, we first compute estimated target size and slope for each trial (each has its own optival)
  # Then with these values we compute partial log-likelihood for "successful" and "failed" trials with the normal CDF, sum to get total loglikelihood
  # optival is sampled to between -1 and 1 because the polynomials work better this way, it is scaled to the provided range before calling simfun
  loglik = function(xx) # xx is a vector with the coefficients for estimate and slope 
   {
    ss = grepl("coefsize",names(xx))
    optival_size  = xx[ss] %*% sapply(trials[,"optival"], function(yy) {yy**(0:(sum(ss)-1))})
    ss = grepl("coefslope",names(xx))
    optival_slope = exp(xx[ss] %*% sapply(trials[,"optival"], function(yy) {yy**(0:(sum(ss)-1))}))
    sum(
         pnorm(intercept_target+optival_slope[!!trials[,"signif"]]*(trials[!!trials[,"signif"],"srsampsize"]-optival_size[!!trials[,"signif"]]),lower.tail=TRUE,log.p=TRUE)
        ,pnorm(intercept_target+optival_slope[ !trials[,"signif"]]*(trials[ !trials[,"signif"],"srsampsize"]-optival_size[ !trials[,"signif"]]),lower.tail=FALSE,log.p=TRUE)
        )
    }
  dloglik = function(xx) # gradient
   {
    # the element-wise gradient calculation by hand shows that it has a few complicated but common components for all arguments (which are different by trial) so let's compute those step by step by trial
    # see paper for the expression of the gradient for slope and size coefficients
    ss_slope = grepl("coefslope",names(xx))
    ss_size = grepl("coefsize",names(xx))
    n_slope = sum(ss_slope)
    n_size = sum(ss_size)
    ordered_vector = xx[c(sort(names(xx)[ss_slope]),sort(names(xx)[ss_size]))] # assuming nobody is going to use a polynomial with x^10 otherwise sorting more complicated
    size_polynomial = colSums(sapply(trials[,"optival"], function(yy) {ordered_vector[paste0("coefsize",0:(n_size-1))]*yy**(0:(n_size-1))}))
    slope_polynomial = colSums(sapply(trials[,"optival"], function(yy) {ordered_vector[paste0("coefslope",0:(n_slope-1))]*yy**(0:(n_slope-1))}))
    optival = trials[,"optival"]
    inner = intercept_target+exp(slope_polynomial)*(trials[,"srsampsize"]-size_polynomial) # 
    signif_sign = c(-1,+1)[1+trials[,"signif"]] # whether upper (+1) or lower (-1) tail of cdf is taken, depending on significant trial or not ; therefore, also sign of dnorm for derivative
    pnorm_inner = sapply(1:nrow(trials),function(yy){pnorm(inner[yy],lower.tail=as.logical(trials[yy,"signif"]))})
    dnorm_inner = dnorm(inner)
    common_slope = dnorm_inner/pnorm_inner*signif_sign*(inner-intercept_target) # common term for gradient terms wrt to slope coeffs
    common_size = dnorm_inner/pnorm_inner*signif_sign*(-1)*exp(slope_polynomial) # common term for gradient terms wrt to size coeffs
    gradient_matrix_distinctive = t(t(optival))[,rep(1,length(xx))]^t(c(0:(n_slope-1),0:(n_size-1)))[rep(1,nrow(trials)),] # picking exponent corresponding to each coefficient
    gradient_matrix_common = cbind(t(t(common_slope))[,rep(1,n_slope)],t(t(common_size))[,rep(1,n_size)])
    gradient = colSums(gradient_matrix_distinctive * gradient_matrix_common)
    names(gradient) = names(ordered_vector)
    gradient = gradient[names(xx)]
    gradient
    }
    # a bit of code to test the gradient #  output = cbind(xx1,xx2) ; for(ii in names(xx1)) { yy = ii ; xx2=xx1 ; xx2[yy]=xx1[yy]+.00000001 ; output[yy,1]=((loglik(xx2)-loglik(xx1))/(xx2[yy]-xx1[yy]))[1] ; output[yy,2]=dloglik(xx1)[yy]} ; print(output) ; print(output[,1]/output[,2])
  
  ddloglik = function(xx) # xx is the same vector, now we compute the Hessian (matrix with second derivatives and derivatives of gradient wrt other parameters)
   {
    # the element-wise Hessian calculation by hand shows that it has a few complicated but common components for all coefficients pairs (and different by trial) so let's compute those step by step by trial
    ss_slope = grepl("coefslope",names(xx))
    ss_size = grepl("coefsize",names(xx))
    n_slope = sum(ss_slope)
    n_size = sum(ss_size)
    ordered_vector = xx[c(sort(names(xx)[ss_slope]),sort(names(xx)[ss_size]))] # assuming nobody is going to use a polynomial with x^10 otherwise sorting more complicated
    size_polynomial = colSums(sapply(trials[,"optival"], function(yy) {ordered_vector[paste0("coefsize",0:(n_size-1))]*yy**(0:(n_size-1))}))
    slope_polynomial = colSums(sapply(trials[,"optival"], function(yy) {ordered_vector[paste0("coefslope",0:(n_slope-1))]*yy**(0:(n_slope-1))}))
    optival = trials[,"optival"]
    inner = intercept_target+exp(slope_polynomial)*(trials[,"srsampsize"]-size_polynomial) # 
    signif_sign = c(-1,+1)[1+trials[,"signif"]] # whether upper (+1) or lower (-1) tail of cdf is taken, depending on significant trial or not ; therefore, also sign of dnorm for derivative
    pnorm_inner = sapply(1:nrow(trials),function(yy){pnorm(inner[yy],lower.tail=as.logical(trials[yy,"signif"]))})
    dnorm_inner = dnorm(inner)
    bigfraction = (-inner*dnorm_inner*pnorm_inner-signif_sign*dnorm_inner^2)/(pnorm_inner^2)
    
    common_size_size = signif_sign*exp(2*slope_polynomial)*bigfraction # common term for the region of the Hessian with two size coeffs
    common_size_slope = -signif_sign*exp(slope_polynomial)*(dnorm_inner/pnorm_inner+(inner-intercept_target)*bigfraction) # common term for the region of the Hessian with one size one slope (formula is symmetrical)
    common_slope_slope = signif_sign*(inner-intercept_target)*(dnorm_inner/pnorm_inner+(inner-intercept_target)*bigfraction) # common term for the region of the Hessian with two slope coeffs
    # term that is not common, is multiplied by the others
    # the distinctive term is, very prosaically, the optival^(sum of the optival exponents corresponding to the coefficients for each pair)
    exponent_matrix = t(c(0:(n_slope-1),0:(n_size-1)))[rep(1,n_slope+n_size),]+t(t(c(0:(n_slope-1),0:(n_size-1)))[rep(1,n_slope+n_size),]) # along the two dimensions of the Hessian, those are the exponents of the shifting term
    hessian_array_distinctive = array(optival,c(nrow(trials),nrow(exponent_matrix),ncol(exponent_matrix)))^array(exponent_matrix,c(1,nrow(exponent_matrix),ncol(exponent_matrix)))[rep(1,nrow(trials)),,]
    hessian_array_common = abind::abind( abind::abind(array(common_slope_slope,c(nrow(trials),n_slope,n_slope)),array(common_size_slope,c(nrow(trials),n_size,n_slope)),along=2)
                                 ,abind::abind(array(common_size_slope,c(nrow(trials),n_slope,n_size)),array(common_size_size,c(nrow(trials),n_size,n_size)),along=2)
                                 ,along=3)
    hessian = apply(hessian_array_common*hessian_array_distinctive,c(2,3),sum)
    dimnames(hessian) = list(rows=names(ordered_vector),cols=names(ordered_vector))
    hessian = hessian[names(xx),names(xx)]
    hessian
    }
    # a bit of code to test the Hessian # output = cbind(xx1,xx2) ; for(ii in names(xx1)) { yy = ii ; xx2=xx1 ; xx2[yy]=xx1[yy]+.00000001 ; output[yy,1]=((dloglik(xx2)-dloglik(xx1))/(xx2[yy]-xx1[yy]))[7] ; output[yy,2]=ddloglik(xx1)[7,yy]} ; print(output) ; print(output[,1]/output[,2]) # repeat with instead of 7 to check another row
    
  optival_rescale = function(xx)        # rescale optimization variable values in xx to the range provided by the user (to be passed to simfun and/or graphs)
   {
    if(optilog) mm = log(optiwin) else mm = optiwin # min and max
    yy = mm[1]+(mm[2]-mm[1])*(xx+1)/2   # (xx+1)/2 scales (-1 to +1) to (0 to 1)
    if(optilog) yy = exp(yy) 
    if(optiround) yy = round(yy)
    unname(yy)
    }
  
  min_NN = NA
  min_optival = NA # best values to date
  # Here is the main loop
  while(  nrow(trials)<nsims                                         # exit if enough simulations made
        & !(nrow(trials)>500 & !is.na(min_NN) & min_NN>impNN)        # exit if after 500+ simulated trials the best size seems to be >impNN
        & !(!is.na(min_NN) & min_NN>1e5 & mean(trials[,"srsampsize"]>sqrt(capNN*.8))>0.1))  # exit in any case if answer seems >100k patients and the simulator is hitting capNN (likely a user error leading to lots of wasted compute)
   {
    if(nrow(trials)>200&&is.null(dim(initiation))&&initiation==FALSE) { trials=trials[-(1:150),] ; initiation = TRUE } # kick out first iterations, deactivate the logical switch
    
    # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # 
    # Running next simulations
    # # # # # # # # # # # # # # # # # # 
    args = list(NN=0,optivar=0,...)  # preparing arguments for vectorized simfun call
    names(args)[2] = optivar         # giving the right name to the argument that is the target of the optimization
    args = lapply(1:length(tar_NN),function(xx){yy=args;yy[[1]]=tar_NN[xx];yy[[2]]=optival_rescale(tar_optival[xx]);yy}) # this is a list of lists, each list has the arguments to be passed to simfun
    
    simulations = unlist(lapply(args,function(xx){do.call(simfun,xx)})) # running the simulations, normally getting a vector of logicals. Note that this calls what you have defined to be simfun when calling adsasi_1d and passes whatever arguments you have specified, except for the optivar parameter which is chosen iteratively by adsasi_1d (values in tar_optival) and sample sizes (values in tar_NN)
    trials = rbind(trials,cbind(srsampsize=sqrt(tar_NN),simulations,optival=tar_optival))  # collating the latest simulations with the ones before
    
    # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # 
    # Fitting model to simulation results
    # # # # # # # # # # # # # # # # # # # # 
    initial_vector = c(c(-5,rep(0,n_slope_coefs-1)),latest_coefs_size)
    names(initial_vector) = c(paste0("coefslope",0:(n_slope_coefs-1)),paste0("coefsize",0:(n_size_coefs-1)))
    initial_value = loglik(initial_vector)                              # log-likelihood for the initial values passed to the minimizer 
    if(initial_value!=Inf&initial_value!=-Inf&!is.na(initial_value))    # if the initial value is defined, the minimizer can get to work
     {
      try(mle <- optim(initial_vector, fn=loglik, gr=dloglik, method="BFGS", hessian = FALSE, control=list(fnscale=-1)), silent = TRUE) # using the default optimizer in R
                                                                                                            # control=list(fnscale=-1)) to maximize instead of minimizing
                                                                                                            # here we use try() in case it fails, gets stuck etc. in which case we keep previous values
      latest_coefs_slope = mle$par[paste0("coefslope",0:(n_slope_coefs-1))]
      latest_coefs_size = mle$par[paste0("coefsize",0:(n_size_coefs-1))]
      
      optival_graph = seq(-1,1,by=.01)                               # 
      slope_by_optival = apply(t(latest_coefs_slope)[rep(1,length(optival_graph)),]*t(t(optival_graph))[,rep(1,n_slope_coefs)]^t(0:(n_slope_coefs-1))[rep(1,length(optival_graph)),],1,sum)^2
      
      covmat = diag(length(mle$par))
      try(covmat <- solve(-ddloglik(c(latest_coefs_slope,latest_coefs_size))), silent=TRUE)                 # using the exact Hessian to get a Rao Cramer covariance matrix, will fail from time to time, mostly when not invertible
      
      # Computing variance of the size polynomial for different values of optival (to get a confidence interval along all values)
      exponent_matrix = t(c(0:(n_slope_coefs-1),0:(n_size_coefs-1)))[rep(1,n_slope_coefs+n_size_coefs),]+t(t(c(0:(n_slope_coefs-1),0:(n_size_coefs-1)))[rep(1,n_slope_coefs+n_size_coefs),]) # along the two dimensions of the Hessian, those are the exponents of the shifting term (powers of v for the covariance of the polynomial)
      
      covar_coefs = array(optival_graph,dim=c(length(optival_graph),1,1))[,rep(1,n_slope_coefs+n_size_coefs),rep(1,n_slope_coefs+n_size_coefs)]^array(exponent_matrix,dim=c(1,nrow(exponent_matrix),ncol(exponent_matrix)))[rep(1,length(optival_graph)),,]
      covar_coefs = abind::abind(covar_coefs,NULL)
      
      covar_array = covar_coefs*array(covmat,dim=c(1,nrow(covmat),ncol(covmat)))[rep(1,length(optival_graph)),,]
      
      size_polynomial_variance = abs(apply(covar_array[,n_slope_coefs+1:n_size_coefs,n_slope_coefs+1:n_size_coefs],1,sum))
      size_se_by_optival = sqrt(size_polynomial_variance)
      size_polynomial_estimate_by_optival = apply(t(latest_coefs_size)[rep(1,length(optival_graph)),]*t(t(optival_graph))[,rep(1,n_size_coefs)]^t(0:(n_size_coefs-1))[rep(1,length(optival_graph)),],1,sum)
      size_confint_lower_by_optival = (size_polynomial_estimate_by_optival-1.96*size_se_by_optival)^2
      size_confint_higher_by_optival = (size_polynomial_estimate_by_optival+1.96*size_se_by_optival)^2
      size_natural_estimate_by_optival = size_polynomial_estimate_by_optival^2
      min_NN = ceiling(min(size_natural_estimate_by_optival))                     # going with an empirical minimum since we have 200 points
      min_optival = optival_graph[which.min(size_natural_estimate_by_optival)[1]] # best optival

      slope_polynomial_variance = abs(apply(covar_array[,1:n_slope_coefs,1:n_slope_coefs],1,sum))
      slope_se_by_optival = sqrt(slope_polynomial_variance)
      slope_polynomial_estimate_by_optival = apply(t(latest_coefs_slope)[rep(1,length(optival_graph)),]*t(t(optival_graph))[,rep(1,n_slope_coefs)]^t(0:(n_slope_coefs-1))[rep(1,length(optival_graph)),],1,sum)
      slope_confint_lower_by_optival = exp(slope_polynomial_estimate_by_optival-1.96*slope_se_by_optival)
      slope_confint_higher_by_optival = exp(slope_polynomial_estimate_by_optival+1.96*slope_se_by_optival)
      slope_natural_estimate_by_optival = exp(slope_polynomial_estimate_by_optival)

      
      batch_size = max(50,round(nrow(trials)/10))                                  # the number of simulations for the new batch is 10% of what's already been done, but at least 50
      if((nrow(trials)+batch_size)>nsims) batch_size = max(2,nsims-nrow(trials))   # max(2,...) is to avoid bugs so one may rarely get nsims+1 runs at the end, depending on nsims
      n_per_optival = rmultinom(1,batch_size,1/size_natural_estimate_by_optival^2) # how many trials per optival value (mostly 0 or 1)
      tar_trials = lapply(1:nrow(n_per_optival),function(xx)                       # this function will pick new locations to try in the {NN,optival} space. We will favor regions near the lowest sample size (+-1SE) and near the best optival (with minimal sample size stimate), but not too harshly (if case we are currently wrong)
       {
        output = NULL
        if(n_per_optival[xx]>0)
         {
          output = round(runif( n_per_optival[xx] 
                               ,min=size_polynomial_estimate_by_optival[xx]-1.96*size_se_by_optival[xx]
                               ,max=size_polynomial_estimate_by_optival[xx]+1.96*size_se_by_optival[xx]
                               )^2)
          output = cbind(output,optival_graph[xx])
          }
        output
        })
      tar_trials = do.call(rbind,tar_trials)      
      
      tar_NN = pmin(pmax(4,tar_trials[,1]),round(runif(nrow(tar_trials),capNN*0.9,capNN)))  # applying limits : minimum N=4, maximum=0.9*capNN to capNN, if hard cap gradient is not estimated properly
      tar_optival = tar_trials[,2]
      } else {                                                                              # if optimizer not launchable, keep the same estimates + some noise, keep target trials (N & optival)
              latest_coefs_size = latest_coefs_size+rnorm(length(latest_coefs_size))
              latest_coefs_slope = latest_coefs_slope+rnorm(length(latest_coefs_slope))
              warning(paste(sep="","Adsasi wrapper log-likelihood fit fail after ",nrow(trials)," simulations."))  # returning a warning, this is usually bad if more than once/twice per run
              }       
        
    if(all(c(trials[,1]^2,tar_NN)==tar_NN[1])) tar_NN = tar_NN*rep(c(.5,1),c(round(length(tar_NN)),length(tar_NN)-round(length(tar_NN))))  # we divide half of values by 2 if all trials have the same size
    if(verbose) { cat("\n") ; print(tar_NN); cat("\n") ; print(se) ; cat("\n") }    # another verbose descriptions of how the run is going
    
    # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # 
    # Plotting
    # # # # # # # # # # # # # # # # # # 
    if(savegraphs!=FALSE) { batch=batch+1 ; pngname = paste(sep="","adsasi ",gsub(":","-",as.character(Sys.time()),fixed=TRUE)," iteration ",paste(rep(0,3-nchar(batch)),collapse=""),batch,".png") ; if(is.character(savegraphs)){pngname=gsub("adsasi",savegraphs,pngname)} ; png(pngname,width=400*3,height=400) } # opening a device to save graphs
    #layout(mat = matrix(c(1,1,1,1,2,3),nrow=2), heights = c(1,1), widths = c(1,1,1), respect =TRUE)    
    par(mfrow=c(1,3),cex=1.35,cex.lab=1,cex.axis=1)  # light aesthetic setup, backup was done before the loop
    on.exit(par(par_bak)) # restoring on exit
    plot(
          trials[,"optival"]                           # optival as abscissae (it will be drawn as -1 to +1 so axis will need to be modified)
         ,trials[,"srsampsize"]^2                      # sample sizes as ordinates 
         ,xlim=c(-1,+1),xaxt="n"                                # optimization space
         ,ylim=list(min_NN*c(.9,2),quantile(trials[,"srsampsize"],c(.01,.99))^2)[[1+is.na(min_NN)]]  # zoomed, not showing everything
         ,main=c(paste(sep="","Status after ",nrow(trials)," simulations"),paste(sep="",c("Current","Final")[1+(nrow(trials)>=nsims)]," minimum ",min_NN," at ",optival_rescale(min_optival))) # writing the latest estimate
         ,xlab=optivar
         ,ylab="Sample size",type="p"
         ,col=paste(sep="",c(rgb(1,.5,0),rgb(0,.6,0.8)),"")[1+trials[,"signif"]],pch=3,cex=.5+.1*(1-trials[,"signif"]) # successes in blue, failures in orange ; so higher ordinates are bluer, lower are oranger
         )
    axis(1,at=c(-1,-.5,0,.5,1),labels=optival_rescale(c(-1,-.5,0,.5,1)))
    if(any(ls()=="size_confint_higher_by_optival"))
     {
      lines(optival_graph,size_natural_estimate_by_optival,lwd=2)
      polygon(c(optival_graph,rev(optival_graph)),c(size_confint_lower_by_optival,rev(size_confint_higher_by_optival)),col="#55555544",border=NA)
      }
    plot(
          tar_optival                                  # optival as abscissae (it will be drawn as -1 to +1 so axis will need to be modified)
         ,tar_NN                                       # sample sizes as ordinates 
         ,xlim=c(-1,+1),xaxt="n"                       # optimization space
         ,ylim=list(min_NN*c(.9,2),quantile(trials[,"srsampsize"],c(.01,.99))^2)[[1+is.na(min_NN)]]  # same as above
         ,main=c(paste(sep="","Next batch of ",length(tar_optival)," simulations"),"(No next batch to be drawn)")[1+(nrow(trials)>=nsims)] # 
         ,xlab=optivar
         ,ylab="Sample size",type="p"
         ,col=c("black",NA)[1+(nrow(trials)>=nsims)] 
         )
    axis(1,at=c(-1,-.5,0,.5,1),labels=optival_rescale(c(-1,-.5,0,.5,1)))
    if(any(ls()=="size_confint_higher_by_optival"))
     {
      lines(optival_graph,size_natural_estimate_by_optival,lwd=2)
      polygon(c(optival_graph,rev(optival_graph)),c(size_confint_lower_by_optival,rev(size_confint_higher_by_optival)),col="#55555544",border=NA)
      }
    
    if(verbose) { print(trials[nrow(trials)-2:0,]) ; cat("\n") ; print(c("estimate_coefs"=latest_coefs_size, "estimate_slopes"=latest_coefs_slope)) }    # verbose descriptions of how the run is going

    if(any(ls()=="slope_confint_higher_by_optival")) { y_window = c(0,quantile(c(slope_natural_estimate_by_optival,slope_confint_lower_by_optival,slope_confint_higher_by_optival),.95)) } else { y_window=c(0,quantile(slope_natural_estimate_by_optival,.95)) }
    plot(                                               # second graph (right panel)
          optival_graph
         ,slope_natural_estimate_by_optival
         ,xlim=c(-1,+1),xaxt="n"                        # optimization space
         ,ylim=y_window                                 # not showing everything
         ,main="Slope",xlab=optivar,ylab="Slope"
         ,type="l",lwd=2,col=2)
    axis(1,at=c(-1,-.5,0,.5,1),labels=optival_rescale(c(-1,-.5,0,.5,1)))
    if(any(ls()=="slope_confint_higher_by_optival")) { polygon(c(optival_graph,rev(optival_graph)),c(slope_confint_lower_by_optival,rev(slope_confint_higher_by_optival)),col="#AA555544",border=NA) }
    if(savegraphs!=FALSE) dev.off()                     # closing device and saving image if graphs are to be saved
    }                                                   # end of loop
  
  if(nrow(trials)<nsims) warning("Halted for inefficiency.")
  if( min_NN > 4*capNN ) { warning("Optimal sample size seems much larger than user-provided size cap, returning +Inf.") ; min_NN=Inf }  # there is little accuracy when trying to extrapolate that far
  output = list(min_NN=min_NN,min_optival=optival_rescale(min_optival))
  if(keepsims) 
   {
    abscissae = optival_rescale(optival_graph)
    trials = cbind(trials, optival_natural = optival_rescale(trials[,"optival"]))
    extra_output = list( trials = trials
                        ,abscissae = abscissae
                        ,latest_coefs_size = latest_coefs_size
                        ,latest_coefs_slope = latest_coefs_slope
                        ,slope_natural_estimate_by_optival = slope_natural_estimate_by_optival
                        ,slope_confint_lower_by_optival = slope_confint_lower_by_optival
                        ,slope_confint_higher_by_optival = slope_confint_higher_by_optival
                        ,size_natural_estimate_by_optival = size_natural_estimate_by_optival
                        ,size_confint_lower_by_optival = size_confint_lower_by_optival
                        ,size_confint_higher_by_optival = size_confint_higher_by_optival)
    output = c(output,extra_output)
    }
  output
  }
