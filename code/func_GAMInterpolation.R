################################################################################################
################################################################################################
## This function estimates the number of outmigrating salmonids using screwtrap capture-recapture
## data and a generalized additive model (GAM) with a penalized spline smooth term.
## The program can use with either continious sample data or subsampled data
## Program Beta Version 2.0 
## By: Adam Duarte & Jim Peterson, Oregon Cooperative Fish and Wildlife Research Unit
## February 1, 2020
##      
## NOTE-->REQUIRES R PACKAGE "mgcv" 
##################################################################################################
## A bit about the uncertainty estimates generated by this function:
##    a) The function estimates confidence intervals for the daily estimates that account for the smoothing term in the GAM. 
##    b) For the group and total abundance estimates, the daily SE estimates as treated as independent and 
##       the total SE for these new intervals is estimated using the delta method. 
##    c) The group and total abundance confidence interval estimates are calculated two ways
##        1) The "Lower_95" and "Upper_95" CI estimates take into account the mean/variance of the point estimate and
##           assumes that the sampling distribution of the estimator is approximately normal. This method does well when
##           the mean estimate is further away from zero. When the mean estimate is near zero, however, the lower CI can sometimes
##           be negative, which is incorrect.
##        2) The "Lower_95_lognorm" and "Upper_95_lognorm" CI estimates also take into account the mean/variance of 
##           the point estimate, but is assumes that the sampling distribution is log-normal (long tail to the right). This
##           method does well when the mean estimate is near zero because it constrains the estimates to be >=0.
##    d) Note that the function sometimes estimates the daily upper CI as infinity. This has to do with the large SEs when there 
##       is no data (i.e., lots of zeroes for catch) and exp(710)=infinite in R.  
##    e) Note that when a group only has one data point the group-level CI's are carried forward from the daily estimates
##################################################################################################

## FUNCTION INPUTS 
### data - name of dataframe containg catch data
### time - name of column containing time index,e.g., day of run
### groups - name of column containing group index,e.g., week. Estimates will be output according to this grouping.
### catch - name of column containing number of fish caught during time interval
### mark - name of column containing number of marked fish released during time interval
### recap - name of column containing number of marked fish recaptured during time interval
### sub.sample - if there are subsampling days in the data (T or F). F is the default.
### P24 - vector of P24 values. NULL is the default.
### PSS - vector of PSS values. NULL is the default.


#### BEGINNING OF FUNCTION
GAM.pSpline.trap.est<-function(data,time ="DAY",catch="CATCH", groups="JWEEK", mark="MRK", recap="RCP",
                               sub.sample=F, sub.sample.data, stime="subDay", pss="PSS", p24="P24", diagnostics=TRUE){
  
  if(is.element("mgcv",installed.packages()[,1])==0){print("ERROR: Package mgcv is not installed");break}
  require(mgcv)
  
  work1<-data
  work<-na.omit(data)  #make sure difference in nrow is number of catch=NA days; see miss below
  
  time=toupper(time);catch=toupper(catch);groups=toupper(groups);mark=toupper(mark);recap=toupper(recap)
  colnam1<-toupper(names(work1))
  colnam<-toupper(names(work))
  
  miss<-nrow(data)-nrow(work)  #number of NA days
  if(miss>0){cat(paste("\nNote",paste(miss," observations excluded due to missing data \n\n",sep=""),sep = " "))}
  
  grp<-work1[,which(colnam1==groups)] #again, grp comes from the complete data frame, not the one filtered to exclude NAs.
  grp0<-work[,which(colnam1==groups)]  #PG edit line. A version of grp that has NAs eliminated. Better wd be to make this grp and prv grp1.
  ch<-work[,which(colnam==catch)]
  mrk<-work[,which(colnam==mark)]
  rcp<-work[,which(colnam==recap)]
  tperiod.1<-work[,which(colnam==time)]
  tperiod=tperiod.1-(min(tperiod.1)-1) # to make the range for time start at 1
  
  if((all(ch==floor(ch)))==0){print("ERROR: Catch data must only contain integer values or NAs");break}
  if((all(mrk==floor(mrk)))==0){print("ERROR: Number marked must only contain integer values");break}
  if((all(rcp==floor(rcp)))==0){print("ERROR: Number recaptured must only contain integer values");break}
  
  bad<-sum(ifelse(rcp>mrk,1,0))
  
  if(bad>0){print(paste("ERROR: Number recaptured is greater than marked for ",paste(bad," observation(s)",sep=""),sep="")); break}
  
  if(sub.sample == T){
    subSampWork<-na.omit(sub.sample.data)
    miss<-nrow(sub.sample.data)-nrow(subSampWork)
    if(miss>0){cat(paste("\nNote",paste(miss," observations excluded from subsample data due to missing data \n\n",sep=""),sep = " "))}
    stime=toupper(stime);pss=toupper(pss);p24=toupper(p24)
    colnamSub<-toupper(names(subSampWork))
    subTime<-subSampWork[,which(colnamSub==stime)]
    subPSS<-subSampWork[,which(colnamSub==pss)]
    subP24<-subSampWork[,which(colnamSub==p24)]
    
    pro<-rep(1,length(work[,1]))
    tempFile1Sub<-data.frame(tperiod,pro)
    colnames(tempFile1Sub)=c("daymatch","pro")
    subDay<-subTime-(min(tperiod.1)-1)  #converting to time index starting at one, I think
    weeklyCal<-(subPSS+1)/(subP24+1)  #I assume plus 1 for same logic as TE?
    tempFile2Sub<-data.frame(subDay,weeklyCal)
    colnames(tempFile2Sub)=c("daymatch","cal")
    tempFile3Sub=merge(tempFile1Sub, tempFile2Sub, by="daymatch",sort=TRUE,all=TRUE)
    PRO<-ifelse(is.na(tempFile3Sub$cal)==TRUE,tempFile3Sub$pro,tempFile3Sub$cal)
    bad<-sum(ifelse(weeklyCal>1,1,0))
    if(bad>0){print(paste("ERROR: Number fish captured in subsample is greater than full sample for ", paste(bad," observation(s)", sep = ""),sep = "")); break}
  } else{PRO=1}
  
  knots<-min(round(max(tperiod)/4),35)
  p=(rcp+1)/(mrk+1)  #chapman version
  est.ch<-round(ch/(p*PRO))  #chapman version.  #estimated abundance, not estimated catch
  tempFile=data.frame(tperiod,ch,p,est.ch)
  colnames(tempFile)=c("tperiod","rawCatch","efficiency","est.Catch")  #again it's abundance not est.Catch!!
  fm<-gam(est.ch~s(tperiod,k=knots,bs="ps"),data=tempFile,family=poisson)
  if (diagnostics==TRUE) {  #PG addition
    gam.check(fm)
  }
  
  outz<-predict.gam(fm,newdata=list(tperiod=1:max(tperiod)),se.fit=TRUE,type="link")
  outz2<-predict.gam(fm,newdata=list(tperiod=1:max(tperiod)),se.fit=TRUE,type="response")
  
  timeStep.mean<-as.vector(outz$fit)
  timeStep.se<-as.vector(outz$se.fit)
  
  timeStep.mean2<-as.vector(outz2$fit)
  timeStep.se2<-as.vector(outz2$se.fit)
  
  rmvn<-function(n,mu,sig) { ## MVN random deviates  #*multivariate normal?
    L<-mroot(sig)
    m<-ncol(L)
    t(mu+L%*%matrix(rnorm(m*n),m,n))
  }
  Vb<-vcov(fm)  #variance-covariance matrix  #dimensions corresponds to number of knots
  
  BUdiff<-rmvn(10000,mu=rep(0,nrow(Vb)),sig=Vb)  #mean is zero, sig=variance covariance matrix. 10k draws from...something? What does BU stand for?
  Cg<-predict(fm,list(tperiod=1:max(tperiod)),type="lpmatrix")  #values of the linear predictor postmltiple by the parameter vector #What does Cg stand for
  simDev<-Cg%*%t(BUdiff) #simulated deviance...Cg times transpose of BUdiff? (matrix multiplication)
  absDev<-abs(sweep(simDev,1,timeStep.se,FUN="/")) #?
  masd<-apply(absDev,2L,max)  #?
  crit<-quantile(masd,probs=0.95,type=8)
  timeStep.down<-fm$family$linkinv(timeStep.mean-(crit*timeStep.se))
  timeStep.up<-fm$family$linkinv(timeStep.mean+(crit*timeStep.se))
  
  tada<-round(data.frame(min(tperiod.1):max(tperiod.1),timeStep.mean2,timeStep.se2,timeStep.down,timeStep.up),2)
  colnames(tada)=c(time,"Estimate","SE","Lower_95","Upper_95")
  # ---- so far that has al been about calculating stats by day
  
  # Starting calc by week[?]
  groupStepData<-data.frame(grp,outz2$fit,outz2$se.fit)
  colnames(groupStepData)= c("group","mean.est","se.est")
  groupStepData$var<-groupStepData$se.est*groupStepData$se.est  #var is std.error squared. Still working with per day values.
  groupStep.mean<-aggregate(groupStepData$mean.est, by=list(Category=groupStepData$group), FUN=sum) #sum daily abundance estimates by week
  groupStep.mean<-groupStep.mean$x #pull out summed weekly abundance values as a vector
  groupStep.var<-aggregate(groupStepData$var, by=list(Category=groupStepData$group), FUN=sum)  #group var as sum of daily var values
  
  n.per.group<-as.vector(table(groupStepData$group))
  n.per.group.names<-as.numeric(names(table(groupStepData$group)))  #I assume this will require that group index be numeric
  thisGroup<-n.per.group.names[which(n.per.group==1)]  #pull out group indexes for groups with only a single daily obs. *This is for the expanded version - should it not have to do with NAs in orig data set? Otherwise only gaps will be start and end - but maybe this is what is meant?
  thisOriginalLine<-which(grp %in% thisGroup) #*I think miscoded, shouldn't be "J_WEEK". Also "data" does not have colnames to all caps.
  thatNewLine<-which(n.per.group==1)
  
  groupStep.se<-sqrt(groupStep.var$x) 
  groupStep.cv<-groupStep.se/groupStep.mean
  C<-exp(1.96*sqrt(log(1+groupStep.cv^2)))  #is this to do with the lognormal version of CI calculation?
  x.lower<-groupStep.mean/C  #I think this is CI version 1
  x.upper<-groupStep.mean*C
  
  groupStep.mean2<-groupStep.mean+groupStep.var$x  #so, why are mean + var added?

  #assume groups have enough individuals in them to approximate Poisson distribution using a normal distribution  #is that individual fish? or individual data points?
  groupStep.down<-groupStep.mean-1.96*sqrt(groupStep.mean2) #.mean2 = mean +var[???]  #CI version 2 [?]
  groupStep.up<-groupStep.mean+1.96*sqrt(groupStep.mean2)
  
  #fix CI for groups that only have 1 observation
  if(length(thisGroup)>1){  #This code replaces the group version with the individual version. so what would the group version look like before this replacement?
    for(kk in 1:length(thisGroup)){
      groupStep.down[thatNewLine[kk]]<-timeStep.down[thisOriginalLine[kk]]
      groupStep.up[thatNewLine[kk]]<-timeStep.up[thisOriginalLine[kk]]
      
      x.lower[thatNewLine[kk]]<-timeStep.down[thisOriginalLine[kk]]
      x.upper[thatNewLine[kk]]<-timeStep.up[thisOriginalLine[kk]]
      
    }
  } else{  #includes length=1 and length=0 (ie, no groups have only one obs). I do think that was a mistake earlier - should be based on data points, not timesteps. 
    groupStep.down[thatNewLine]<-timeStep.down[thisOriginalLine]
    groupStep.up[thatNewLine]<-timeStep.up[thisOriginalLine]
    
    x.lower[thatNewLine]<-timeStep.down[thisOriginalLine]
    x.upper[thatNewLine]<-timeStep.up[thisOriginalLine]
  }
  
  tada2<-round(data.frame(min(grp):max(grp),groupStep.mean,groupStep.se,groupStep.down,groupStep.up,x.lower,x.upper),2) #This code assumes group index is continusous. I think it wd be more robust to just select distinct
  colnames(tada2)=c(groups,"Estimate","SE","Lower_95","Upper_95","Lower_95_lognorm","Upper_95_lognorm")
  
  total.mean<-sum(outz2$fit)
  total.se<-sqrt(sum(outz2$se.fit^2)) 
  total.var<-sum(outz2$se.fit^2) 
  total.mean2<-total.mean+total.var #??
  total.cv<-total.se/total.mean
  CTotal<-exp(1.96*sqrt(log(1+total.cv^2)))
  total.lower<-total.mean/CTotal
  total.upper<-total.mean*CTotal

  #assume there is enough individuals total to approximate Poisson distribution using a normal distribution
  total.low<-total.mean-1.96*sqrt(total.mean2)
  total.high<-total.mean+1.96*sqrt(total.mean2)

  tada3<-round(data.frame(total.mean,total.se,total.low,total.high,total.lower,total.upper),2)
  colnames(tada3)=c("Estimate","SE","Lower_95","Upper_95","Lower_95_lognorm","Upper_95_lognorm")
  
  output<-list(tada,tada2,tada3)
  names(output)<-c("Time_specific_estimates","Group_specific_estimates","Total_population_estimate")
  return(output)
}
##### END OF FUNCTION

## FUNCTION INPUTS 
### data - name of dataframe containg catch data
### time - name of column containing time index,e.g., day of run
### groups - name of column containing group index,e.g., week. Estimates will be output according to this grouping.
### catch - name of column containing number of fish caught during time interval
### mark - name of column containing number of marked fish released during time interval
### recap - name of column containing number of marked fish recaptured during time interval
### sub.sample - if there are subsampling days in the data (T or F). F is the default.
### P24 - vector of P24 values. NULL is the default.
### PSS - vector of PSS values. NULL is the default.