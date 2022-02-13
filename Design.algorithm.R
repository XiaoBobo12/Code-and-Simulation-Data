library(markovchain)
##################################
## simulation procedures
##################################

hazard.function.control<-function(t,k,h0)
  ## the hazard function of the control group
  ## t:the survival time
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
{
  return(k*h0*(h0*t)^(k-1))
}

hazard.function.treatment<-function(t,k,h0,t1,t2,theta,a,b)
  ## the hazard function of the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t:the survival time
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  if(t1==t2)
  {
    h0t=k*h0*(h0*t)^(k-1)
    h1t=theta*h0t
    return((t<t1)*h0t+(t>=t1)*h1t)
  }
  else
  {
    h0t=k*h0*(h0*t)^(k-1)
    t.std=(t-t1)/(t2-t1)
    lt=0*(t<=t1)+pbeta(t.std,a,b)*(t1<t && t<=t2)+(t>t2)
    h1t=(1-lt+theta*lt)*h0t
  }
  return(h1t)
}

cum.hazard.function.treatment<-function(t,k,h0,t1,t2,theta,a,b)
  ## the cumulative hazard function of the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  cum.hazard.treatment.unit<-integrate(f=hazard.function.treatment,lower=0,upper=t,k,h0,t1,t2,theta,a,b)$value
  return(cum.hazard.treatment.unit)
}

cum.hazard.function.treatment.ref<-function(t,k,h0,t1,t2,theta,a,b,Cum.Hazd)
  ## the standardialized cumulative hazard function of the treatment group, which is used to produce the survival time
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  return(cum.hazard.function.treatment(t,k,h0,t1,t2,theta,a,b)-Cum.Hazd)
}

survival.function.control<-function(t,k,h0)
  ## the survival function of the control group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
{
  return(exp(-(h0*t)^k))
}

survival.function.treatment<-function(t,k,h0,t1,t2,theta,a,b)
  ## the survival function of the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  if(length(t)==1)
  {
    survival.proportion=exp(-(cum.hazard.function.treatment(t,k,h0,t1,t2,theta,a,b)))
    return(survival.proportion)
  }
  else
  {
    n=length(t)
    survival.proportion<-rep(NA,n)
    for(i in 1:n)
    {
      survival.proportion[i]=exp(-(cum.hazard.function.treatment(t[i],k,h0,t1,t2,theta,a,b)))
    }
    return(survival.proportion)
  }
}

survival.time.control<-function(n,k,h0)
  ## -n:number of observations. If length(n) > 1, the length is taken to be the number required
  ## -rate0:rates before the time of change
  ## -k:the weibull shape parameter
{
  U=runif(n)
  radnm=(-log(U))^(1/k)/h0
  return(radnm)
}

survival.time.treatment<-function(n,k,h0,t1,t2,theta,a,b)
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  if(theta==1)
  {
    U=runif(n)
    survival.time=(-log(U))^(1/k)/h0
  }
  else
  {
    U=runif(n)
    Cum.Hazd=-log(U)
    survival.time<-rep(NA,n)
    for(i in 1:n)
    {
      survival.time[i]=uniroot(f=cum.hazard.function.treatment.ref,interval=c(0,10000),k,h0,t1,t2,theta,a,b,Cum.Hazd[i])$root
    }
  }
  return(list(survival.time=survival.time,survival.proportion=U))
}


rpoipross<-function(n,rate)
  ## produce the entry time points accroding to the Poisson process
  ## -n:the number of time points
  ## -rate:intensity
{
  if(length(n)!=1)
    stop("n is not one parameter.")
  if(length(rate)!=1)
    stop("rate is not one parameter.")
  if(is.integer(n))
  {
    stop("n is not integer!")
  }
  if(n<=0)
  {
    stop("n is not bigger than 0")
  }
  if(!is.double(rate))
  {
    stop("rate is not a decimal or integer")
  }
  if(rate<=0)
  {
    stop("rate is not bigger than 0")
  }
  t=0
  S=rep(0,n)
  for(i in c(1:n))
  {
    U<-runif(1)
    t=t-log(U)/rate
    S[i]=t
  }
  return (S)
}

dataset.produce<-function(nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1,t2,theta,a,b)
  ## produce the dataset used for the hypothesis in the Monte-Carlo simulation-procedure
  ## nc,nt:the sample sizes of the control and treatment groups
  ## accrual.type:select the way of how patients enter ths study
  ## c(0):patients enter the study according to the poisson process
  ## c(1,a):patients enter the study according to the distribution of F(x)=(x/A)^a
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## split.num: the splitted number of an unit period for an Markov chain
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  datasetc<-matrix(nrow=nc,ncol=7)
  datasett<-matrix(nrow=nt,ncol=7)
  datasetc[,1]=rep(1,nc)## "1" represents the control group 
  datasett[,1]=rep(0,nt)## "0" represents the treatment group 
  h.loss.control=-log(1-loss.control)## the hazard constant of loss to follow-up in the control group
  h.loss.treatment=-log(1-loss.treatment)## the hazard constant of loss to follow-up in the treatment group
  ##produce the initial survival time of two groups
  datasetc[,2]=survival.time.control(nc,k,h0)
  datasett[,2]=survival.time.treatment(nt,k,h0,t1,t2,theta,a,b)$survival.time
  ##produce the recruitment time of two groups
  if(accrual.type[1]==0)
  {
    ##produce the entry time of control group and treatment group by poisson process
    entry=rpoipross(nc+nt,(nc+nt)/A)
    entry.total.order=c(1:(nc+nt))## original order number of two groups
    entry.control.order=sample(entry.total.order,nc)## select the order number corresponding to the entry time of the control group randomly
    entry.treatment.order=entry.total.order[-match(entry.control.order,entry.total.order)]## the rest is the order number of the treatment group
    entry.treatment.order=entry.treatment.order[sample(1:nt)]## permutate the order number of the treatment group
    datasetc[,3]=entry[entry.control.order]## get the entry number of the control group
    datasett[,3]=entry[entry.treatment.order]## get the entry number of the treatment group
  }
  else
  {
    a=accrual.type[2]
    datasetc[,3]=A*(runif(nc,0,1))^(1/a)
    datasett[,3]=A*(runif(nt,0,1))^(1/a)
  }
  ## produce the random censoring time
  datasetc[,4]=rexp(nc,h.loss.control)
  datasett[,4]=rexp(nt,h.loss.treatment)
  ## if the survial time is longer than the  censor time,the subject is censored;uncersored vice versa
  datasetc[,5]= datasetc[,2]<datasetc[,4]
  datasett[,5]= datasett[,2]<datasett[,4]
  ## obtain the observed time since the recruitment
  datasetc[,6]=apply(cbind(datasetc[,2],datasetc[,4]),1,min)
  datasett[,6]=apply(cbind(datasett[,2],datasett[,4]),1,min)
  ## obtain the observed time since the beginning of the clinical trial
  datasetc[,7]=datasetc[,6]+datasetc[,3]
  datasett[,7]=datasett[,6]+datasett[,3]
  dataset.whole=matrix(nrow=nc+nt,ncol=8)## the whole dataset of the control arm and the treatment arm 
  dataset.cb0=rbind(datasett,datasetc)## combine two datasets
  dataset.cb=dataset.cb0[order(dataset.cb0[,7]),]## sort the matrix by the column,which represents "the end time of actual visit"
  dataset.whole[,1]=seq(1,nc+nt)## the order number of two groups
  dataset.whole[,2:8]=dataset.cb##the data need to store
  colnames(dataset.whole)=c("order.number","group","survival.time","entry.time","censorring.time","censor.or.death","actual.observed.time.since.recruitment","actual.observed.time.since.begining")
  return(dataset.whole)
}

log.rank.test<-function(dataset.input)
  #This function is used to calculate the value of the log-rank test statistic
{
  dataset.analysis<-dataset.input[,c(2,6,7)]## group, censor.or.death, actual.observed.time.since.recruitment
  dataset.analysis<-dataset.analysis[order(dataset.analysis[,3]),]
  dataset.uncensored<-dataset.analysis[dataset.analysis[,2]==1,]## remove the censored patients
  number.uncensored<-length(dataset.uncensored[,1])## get the number of uncensored subject
  nc<-sum(dataset.analysis[,1]==1)
  ne<-sum(dataset.analysis[,1]==0)
  ncvector=rep(0,nc+ne+1)
  pcvector=rep(0,nc+ne+1)
  nevector=rep(0,nc+ne+1)
  pevector=rep(0,nc+ne+1)
  ncvector[1]=nc
  nevector[1]=ne
  pcvector[1]=nc/(nc+ne)
  pevector[1]=ne/(nc+ne)
  Sw=matrix(0,nrow=number.uncensored,ncol=3)
  ## the data used to calculate log-rank statistics 
  ## the first row represents the order number of the death after the onset of effect
  ## the second row represents the  numerator of the log-rank statistics
  ## the third row represents the square of the denumerator of the log-rank statistics
  iterator=0## an iterator used to calculate the log-rank statistics,which represents the number of deaths used in the calculation of log-rank test statistics
  for(i in 1:(nc+ne+1))
  {
    if(dataset.analysis[i,1]==1)## the number of the patiens in control group minus 1,while the other remains unchanged
    {
      ncvector[i+1]=ncvector[i]-1
      nevector[i+1]=nevector[i]
    }
    if(dataset.analysis[i,1]==0)## the number of the patiens in treatment group minus 1,while the other remains unchanged
    {
      ncvector[i+1]=ncvector[i]
      nevector[i+1]=nevector[i]-1
    }
    if(ncvector[i+1]==0 & nevector[i+1]==0)## once the number of either group is 0,stop the calculation of log-rank statistics
      break
    pcvector[i+1]=ncvector[i+1]/(ncvector[i+1]+nevector[i+1])
    pevector[i+1]=nevector[i+1]/(ncvector[i+1]+nevector[i+1])
    if(dataset.analysis[i,2]==1)
    {
      ## uncesored
      iterator=iterator+1##iterator add 1
      Sw[iterator,1]=iterator
      if(iterator!=1)
      {
        Sw[iterator,2]=Sw[iterator-1,2]+(dataset.analysis[i,1]- pcvector[i])
        Sw[iterator,3]=Sw[iterator-1,3]+pcvector[i]*pevector[i]
      }
      else
      {
        Sw[iterator,2]=(dataset.analysis[i,1]- pevector[i])
        Sw[iterator,3]=pcvector[i]*pevector[i]
      }
    }
  }
  ##get the value of the log-rank statistics corresponding to interim analysis and final analysis
  Sw_used=Sw[iterator,]
  Sw_value<-Sw_used[2]/sqrt(Sw_used[3])## calculate the value of the piecewise weitghted log-rank statistics
  return(list(Sw_value=Sw_value,death.all=number.uncensored,num.recru=(nc+ne),Sw.ustd=Sw_used[2],SW.Var.H0=Sw_used[3]))
}

log.rank.Monte.Carlo.simulation<-function(order,nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1.true,t2.true,theta,a,b)
  ## Simulate a fixed sample trial and calculate the log-rank test statistics at the given calendar driven (parallel version)
  ## order: the order number of simulation
  ## nc,nt:the sample sizes of the control and treatment groups
  ## accrual.type:select the way of how patients enter ths study
  ## c(0):patients enter the study according to the poisson process
  ## c(1,a):patients enter the study according to the distribution of F(x)=(x/A)^a
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## noncomp.control: the proportion of the non-compliance subjects in the control group in an unit period
  ## noncomp.treatment: the proportion of the non-compliance subjects in the treatment group in an unit period
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time(used in the MRET test)
  ## t1.true,t2.true:the minimum and maximum delay time (used in the generation of the dataset)
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  analysis.time=tau
  dataset.whole<-dataset.produce(nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1.true,t2.true,theta,a,b)
  dataset.analysis<-dataset.whole[dataset.whole[,4]<=analysis.time,]
  dataset.analysis[dataset.analysis[,8]>analysis.time,6]=0
  dataset.analysis[dataset.analysis[,8]>analysis.time,7]=analysis.time-dataset.analysis[dataset.analysis[,8]>analysis.time,4]
  dataset.analysis[dataset.analysis[,8]>analysis.time,8]=analysis.time
  test.statistic=log.rank.test(dataset.input=dataset.analysis)
  Sw_value=test.statistic$Sw_value
  Sw.ustd=test.statistic$Sw.ustd
  SW.Var.H0=test.statistic$SW.Var.H0
  Events.number.total=test.statistic$death.all
  return(c(order=order,Sw_value=Sw_value,Sw.ustd=Sw.ustd,SW.Var.H0=SW.Var.H0,Events.number.total=Events.number.total))
}

Maximin.efficiency.robust.test<-function(dataset.input,t.low,t.upper)
  ## This function is used to calculate the value of the Maximin efficiency robust test
  ## dataset.input: the dataset used to survival analysis
  ## t.low: the minimum delay time
  ## t.upper: the maximum delay time
{
  dataset.analysis<-dataset.input[,c(2,6,7)]## group, censor.or.death, actual.observed.time.since.recruitment
  dataset.analysis<-dataset.analysis[order(dataset.analysis[,3]),]
  number.uncensored<-sum(dataset.analysis[,2]==1)## get the number of uncensored subject
  nc<-sum(dataset.analysis[,1]==1)
  ne<-sum(dataset.analysis[,1]==0)
  nc.tir<-nc
  ne.tir<-ne
  Calc.matrix<-matrix(nrow=nc+ne,ncol=11)
  colnames(Calc.matrix)<-c("pc","pe","pc.pe","psi","censor.or.not","t1_t2","greater.t2","Group","weight","Ustd.unit","Variance.unit")
  Calc.matrix[1,1]=nc.tir/(nc.tir+ne.tir)
  Calc.matrix[1,2]=ne.tir/(nc.tir+ne.tir)
  Calc.matrix[1,3]=nc.tir*ne.tir/(nc.tir+ne.tir)^2
  Calc.matrix[1,4]=Calc.matrix[1,3]*dataset.analysis[1,2]
  Calc.matrix[1,5]=dataset.analysis[1,2]
  Time.obs=dataset.analysis[1,3]
  Calc.matrix[1,6]=Time.obs>=t.low & Time.obs<=t.upper
  Calc.matrix[1,7]=Time.obs>t.upper
  Calc.matrix[1,8]=dataset.analysis[1,1]
  nc.tir=nc.tir-dataset.analysis[1,1]
  ne.tir=ne.tir-(1-dataset.analysis[1,1])
  for(i in 2:(nc+ne))
  {
    Calc.matrix[i,1]=nc.tir/(nc.tir+ne.tir)
    Calc.matrix[i,2]=ne.tir/(nc.tir+ne.tir)
    Calc.matrix[i,3]=nc.tir*ne.tir/(nc.tir+ne.tir)^2
    Calc.matrix[i,4]=Calc.matrix[i,3]*dataset.analysis[i,2]+Calc.matrix[i-1,4]
    Calc.matrix[i,5]=dataset.analysis[i,2]
    Time.obs=dataset.analysis[i,3]
    Calc.matrix[i,6]=Time.obs>=t.low & Time.obs<=t.upper
    Calc.matrix[i,7]=Time.obs>t.upper
    Calc.matrix[i,8]=dataset.analysis[i,1]
    nc.tir=nc.tir-dataset.analysis[i,1]
    ne.tir=ne.tir-(1-dataset.analysis[i,1])
  }
  Calc.matrix[,4]=Calc.matrix[,4]/(nc+ne)
  psi.tau<-Calc.matrix[nc+ne,4]
  if(sum(Calc.matrix[,6])==0)
  {
    psi.t1=0
    psi.t2=0
  }
  else
  {
    psi.t1<-ifelse(min(which(Calc.matrix[,6]==1))==1,0,Calc.matrix[min(which(Calc.matrix[,6]==1))-1,4])
    psi.t2<-Calc.matrix[max(which(Calc.matrix[,6]==1)),4]
  }
  Calc.matrix[,9]=sqrt((psi.tau-psi.t1)/(psi.tau-Calc.matrix[,4]*Calc.matrix[,6]))*Calc.matrix[,6]+2*sqrt((psi.tau-psi.t1)/(psi.tau-psi.t2))*Calc.matrix[,7]
  Calc.matrix[,10]=Calc.matrix[,9]*(Calc.matrix[,8]-Calc.matrix[,1])*Calc.matrix[,5]
  Calc.matrix[,11]=Calc.matrix[,9]^2*Calc.matrix[,3]*Calc.matrix[,5]
  Value.unstd<-sum(Calc.matrix[,10])
  var.H0<-sum(Calc.matrix[,11])
  Sw_value<-Value.unstd/sqrt(var.H0)
  return(list(Sw_value=Sw_value,death.all=number.uncensored,num.recru=(nc+ne),Sw.ustd=Value.unstd,SW.Var.H0=var.H0,Calc.matrix=Calc.matrix,psi.t2=psi.t2,psi.t1=psi.t1,psi.tau=psi.tau))
}

MERT.Monte.Carlo.simulation<-function(order,nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1,t2,t1.true,t2.true,theta,a,b)
  ## Simulate a fixed sample trial and calculate the maximin efficiency robust test statistics at the given calendar driven (parallel version)
  ## order: the order number of simulation
  ## nc,nt:the sample sizes of the control and treatment groups
  ## accrual.type:select the way of how patients enter ths study
  ## c(0):patients enter the study according to the poisson process
  ## c(1,a):patients enter the study according to the distribution of F(x)=(x/A)^a
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## noncomp.control: the proportion of the non-compliance subjects in the control group in an unit period
  ## noncomp.treatment: the proportion of the non-compliance subjects in the treatment group in an unit period
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time(used in the MRET test)
  ## t1.true,t2.true:the minimum and maximum delay time (used in the generation of the dataset)
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  analysis.time=tau
  dataset.whole<-dataset.produce(nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1.true,t2.true,theta,a,b)
  dataset.analysis<-dataset.whole[dataset.whole[,4]<=analysis.time,]
  dataset.analysis[dataset.analysis[,8]>analysis.time,6]=0
  dataset.analysis[dataset.analysis[,8]>analysis.time,7]=analysis.time-dataset.analysis[dataset.analysis[,8]>analysis.time,4]
  dataset.analysis[dataset.analysis[,8]>analysis.time,8]=analysis.time
  test.statistic=Maximin.efficiency.robust.test(dataset.input=dataset.analysis,t.low=t1,t.upper=t2)
  Sw_value=test.statistic$Sw_value
  Sw.ustd=test.statistic$Sw.ustd
  SW.Var.H0=test.statistic$SW.Var.H0
  Events.number.total=test.statistic$death.all
  return(c(order=order,Sw_value=Sw_value,Sw.ustd=Sw.ustd,SW.Var.H0=SW.Var.H0,Events.number.total=Events.number.total))
}
  
fixed.sample.size.trial.Monte.Carlo.simulation<-function(nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1,t2,t1.true,t2.true,theta,a,b,sides,alpha,test.type,simnum)
  # this function is used to simulate the fixed sample trial
  ## nc,nt:the sample sizes of the control and treatment groups
  ## accrual.type:select the way of how patients enter ths study
  ## c(0):patients enter the study according to the poisson process
  ## c(1,a):patients enter the study according to the distribution of F(x)=(x/A)^a
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## noncomp.control: the proportion of the non-compliance subjects in the control group in an unit period
  ## noncomp.treatment: the proportion of the non-compliance subjects in the treatment group in an unit period
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time(used in the MRET test;default to t1.true and t2.true for the log-rank test)
  ## t1.true,t2.true:the minimum and maximum delay time (used in the generation of the dataset)
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time  
  ## alpha: the significance level
  ## sides: 1:right-sided test;2:two-sided test
  ## test.type:1:log-rank test;2:Maximin efficiency robustness test
  ## simnum: the number of repetitions in the Monte-Carlo simulation procedure
{
  if(sides==1)## the right side
  {
    lowerbound=-Inf## get the lower boundary of the group sequential design
    upperbound=qnorm(1-alpha)## get the upper boundary of the group sequential design
  }
  else if(sides==2)## the two sides
  {
    upperbound<-qnorm(1-alpha/2) ## get the upper boundary of the group sequential design
    lowerbound<--qnorm(1-alpha/2) ## get the lower boundary of the group sequential design
  }
  analysis.time=tau
  test.statistics.simulation=rep(0,simnum)
  Decision.simulation=rep(0,simnum)
  Sw.ustd.simulation=rep(0,simnum)
  SW.Var.H0.simulation=rep(0,simnum)
  Events.number.total.simulation=rep(0,simnum)
  if(test.type==1)
  {
    ## log-rank test
    for(l in 1:simnum)
    {
      print(paste(l,"/",simnum))
      dataset.whole<-dataset.produce(nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1.true,t2.true,theta,a,b)
      dataset.analysis<-dataset.whole[dataset.whole[,4]<=analysis.time,]
      dataset.analysis[dataset.analysis[,8]>analysis.time,6]=0
      dataset.analysis[dataset.analysis[,8]>analysis.time,7]=analysis.time-dataset.analysis[dataset.analysis[,8]>analysis.time,4]
      dataset.analysis[dataset.analysis[,8]>analysis.time,8]=analysis.time
      test.statistic=log.rank.test(dataset.analysis)
      Sw_value=test.statistic$Sw_value
      test.statistics.simulation[l]=Sw_value
      Sw.ustd.simulation[l]=test.statistic$Sw.ustd
      SW.Var.H0.simulation[l]=test.statistic$SW.Var.H0
      Decision.simulation[l]=!(Sw_value<=upperbound&lowerbound<=Sw_value)## the rejection situation for respective situation,1 represents rejection,0 represents acceptance
      Events.number.total.simulation[l]=test.statistic$death.all
    }
  }
  else
  {
    ## MRET test
    for(l in 1:simnum)
    {
      print(paste(l,"/",simnum))
      dataset.whole<-dataset.produce(nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1.true,t2.true,theta,a,b)
      dataset.analysis<-dataset.whole[dataset.whole[,4]<=analysis.time,]
      dataset.analysis[dataset.analysis[,8]>analysis.time,6]=0
      dataset.analysis[dataset.analysis[,8]>analysis.time,7]=analysis.time-dataset.analysis[dataset.analysis[,8]>analysis.time,4]
      dataset.analysis[dataset.analysis[,8]>analysis.time,8]=analysis.time
      test.statistic=Maximin.efficiency.robust.test(dataset.input=dataset.analysis,t.low=t1,t.upper=t2)
      Sw_value=test.statistic$Sw_value
      test.statistics.simulation[l]=Sw_value
      Sw.ustd.simulation[l]=test.statistic$Sw.ustd
      SW.Var.H0.simulation[l]=test.statistic$SW.Var.H0
      Decision.simulation[l]=!(Sw_value<=upperbound&lowerbound<=Sw_value)## the rejection situation for respective situation,1 represents rejection,0 represents acceptance
      Events.number.total.simulation[l]=test.statistic$death.all
    }
  }
  Power=mean(Decision.simulation)## the emprical power
  analysis.time.pred=analysis.time
  Events.pred=mean(Events.number.total.simulation)
  Mean.H1.std=mean(test.statistics.simulation)
  Var.H1.std=var(test.statistics.simulation)
  Mean.H1.ustd=mean(Sw.ustd.simulation)
  Var.H1.ustd=var(Sw.ustd.simulation)
  Var.H0.ustd=mean(SW.Var.H0.simulation)
  Var.Var.H0.ustd=var(SW.Var.H0.simulation)
  Empirical.Distribution=list(Mean.H1.std=Mean.H1.std,Var.H1.std=Var.H1.std,Mean.H1.ustd=Mean.H1.ustd,Var.H1.ustd=Var.H1.ustd,Var.H0.ustd=Var.H0.ustd,Var.Var.H0.ustd=Var.Var.H0.ustd,log.rank.ustd=Sw.ustd.simulation,Var.H0.simulation=SW.Var.H0.simulation)
  Initial.data=list(test.statistics.simulation=test.statistics.simulation,Sw.ustd.simulation=Sw.ustd.simulation,SW.Var.H0.simulation=SW.Var.H0.simulation)
  Parameters.list=list(nc=nc,nt=nt,accrual.type=accrual.type,A=A,tau=tau,loss.control=loss.control,loss.treatment=loss.treatment,k=k,h0=h0,t1=t1,t2=t2,t1.true=t1.true,t2.true=t2.true,theta=theta,a=a,b=b,sides=sides,alpha=alpha,test.type=test.type,simnum=simnum)
  return(list(Power=Power,analysis.time.pred=analysis.time.pred,Events.pred=Events.pred,Parameters.list=Parameters.list,Empirical.Distribution=Empirical.Distribution,Initial.data=Initial.data))
}

library(parallel)
fixed.sample.size.trial.Monte.Carlo.simulation.parallel<-function(nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1,t2,t1.true,t2.true,theta,a,b,sides,alpha,test.type,simnum)
  # this function is used to simulate the fixed sample trial(parallel version)
  ## nc,nt:the sample sizes of the control and treatment groups
  ## accrual.type:select the way of how patients enter ths study
  ## c(0):patients enter the study according to the poisson process
  ## c(1,a):patients enter the study according to the distribution of F(x)=(x/A)^a
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## noncomp.control: the proportion of the non-compliance subjects in the control group in an unit period
  ## noncomp.treatment: the proportion of the non-compliance subjects in the treatment group in an unit period
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time(used in the MRET test;default to t1.true and t2.true for the log-rank test)
  ## t1.true,t2.true:the minimum and maximum delay time (used in the generation of the dataset)
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time  
  ## alpha: the significance level
  ## sides: 1:right-sided test;2:two-sided test
  ## test.type:1:log-rank test;2:Maximin efficiency robustness test
  ## simnum: the number of repetitions in the Monte-Carlo simulation procedure
{
  if(sides==1)## the right side
  {
    lowerbound=-Inf## get the lower boundary of the group sequential design
    upperbound=qnorm(1-alpha)## get the upper boundary of the group sequential design
  }
  else if(sides==2)## the two sides
  {
    upperbound<-qnorm(1-alpha/2) ## get the upper boundary of the group sequential design
    lowerbound<--qnorm(1-alpha/2) ## get the lower boundary of the group sequential design
  }
  analysis.time=tau
  core.num=detectCores()-1
  varlist=c("log.rank.test","Maximin.efficiency.robust.test","hazard.function.control","hazard.function.treatment","cum.hazard.function.treatment","cum.hazard.function.treatment.ref","survival.function.control","survival.function.treatment","survival.time.control","survival.time.treatment","rpoipross","dataset.produce","log.rank.Monte.Carlo.simulation","MERT.Monte.Carlo.simulation")
  if(test.type==1)
  {
    ## log-rank test
    cl<-makeCluster(getOption("cl.cores",core.num))
    clusterExport(cl,varlist=varlist,envir=environment())
    simul.res.list=parLapply(cl,seq(1,simnum,1),log.rank.Monte.Carlo.simulation,nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1.true,t2.true,theta,a,b)
    stopCluster(cl)
    Simul.res <- do.call(rbind,simul.res.list)
    test.statistics.simulation<-Simul.res[,2]
    Decision.simulation=!(Simul.res[,2]<=upperbound&lowerbound<=Simul.res[,2])
    Sw.ustd.simulation=Simul.res[,3]
    SW.Var.H0.simulation=Simul.res[,4]
    Events.number.total.simulation=Simul.res[,5]
  }
  else
  {
    ## MRET test
    cl<-makeCluster(getOption("cl.cores",core.num))
    clusterExport(cl,varlist=varlist,envir=environment())
    simul.res.list=parLapply(cl,seq(1,simnum,1),MERT.Monte.Carlo.simulation,nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1,t2,t1.true,t2.true,theta,a,b)
    stopCluster(cl)
    Simul.res <- do.call(rbind,simul.res.list)
    test.statistics.simulation<-Simul.res[,2]
    Decision.simulation=!(Simul.res[,2]<=upperbound&lowerbound<=Simul.res[,2])
    Sw.ustd.simulation=Simul.res[,3]
    SW.Var.H0.simulation=Simul.res[,4]
    Events.number.total.simulation=Simul.res[,5]
  }
  Power=mean(Decision.simulation)## the emprical power
  analysis.time.pred=analysis.time
  Events.pred=mean(Events.number.total.simulation)
  Mean.H1.std=mean(test.statistics.simulation)
  Var.H1.std=var(test.statistics.simulation)
  Mean.H1.ustd=mean(Sw.ustd.simulation)
  Var.H1.ustd=var(Sw.ustd.simulation)
  Var.H0.ustd=mean(SW.Var.H0.simulation)
  Var.Var.H0.ustd=var(SW.Var.H0.simulation)
  Empirical.Distribution=list(Mean.H1.std=Mean.H1.std,Var.H1.std=Var.H1.std,Mean.H1.ustd=Mean.H1.ustd,Var.H1.ustd=Var.H1.ustd,Var.H0.ustd=Var.H0.ustd,Var.Var.H0.ustd=Var.Var.H0.ustd,log.rank.ustd=Sw.ustd.simulation,Var.H0.simulation=SW.Var.H0.simulation)
  Initial.data=list(test.statistics.simulation=test.statistics.simulation,Sw.ustd.simulation=Sw.ustd.simulation,SW.Var.H0.simulation=SW.Var.H0.simulation)
  Parameters.list=list(nc=nc,nt=nt,accrual.type=accrual.type,A=A,tau=tau,loss.control=loss.control,loss.treatment=loss.treatment,k=k,h0=h0,t1=t1,t2=t2,t1.true=t1.true,t2.true=t2.true,theta=theta,a=a,b=b,sides=sides,alpha=alpha,test.type=test.type,simnum=simnum)
  return(list(Power=Power,analysis.time.pred=analysis.time.pred,Events.pred=Events.pred,Parameters.list=Parameters.list,Empirical.Distribution=Empirical.Distribution,Initial.data=Initial.data))
}


marcov.survival.states.chain.simple<-function(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t2,theta,a,b)
  ## the function used to simulate the state proportions at different time points
  ## Note: this marcov model is not applicable to the switch treatments under the delayed treatment effect pattern
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## split.interval.num: the number of the further splitted small intervals for the interval which needs further
  ## interval.period: the duration of the interval which will be further splitted
  ## noncomp.control: the proportion of the non-compliance subjects in the control group in an unit period
  ## noncomp.treatment: the proportion of the non-compliance subjects in the treatment group in an unit period
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  split.num<-1/interval.period*split.interval.num## the splitted number of the unit time interval for an Markov chain
  ##this alogorithm restricts that split.interval.num is exactly divided by interval.period
  interval.num=tau/interval.period## the total number of the intervals which are initially splitted
  unit.interval.period=interval.period/split.interval.num## the period of the smallest splitted interval
  num<-split.num*tau## the total number of the splitted intervals
  loss.control.split<-1-(1-loss.control)^(1/split.num)
  loss.treatment.split<-1-(1-loss.treatment)^(1/split.num)
  Control.matrix<-matrix(0,nrow=(num+1),ncol=4)
  colnames(Control.matrix)<-c("Time","Loss","Event","Surv.cont")
  Control.matrix[1,4]=1
  Treatment.matrix<-matrix(0,nrow=(num+1),ncol=4)
  colnames(Treatment.matrix)<-c("Time","Loss","Event","Surv.Trt")
  Treatment.matrix[1,4]=1
  Trans.matrix.list.control<-list()
  Trans.matrix.list.treatment<-list()
  for(l in 1:interval.num)
  {
    t.low=(l-1)*interval.period
    t.high=l*interval.period
    cum.hazard.control<-integrate(f=hazard.function.control,lower=t.low,upper=t.high,k,h0)$value
    Prop.death.control.split<-1-exp(-cum.hazard.control/split.interval.num)
    cum.hazard.treatment.split<-integrate(f=hazard.function.treatment,lower=t.low,upper=t.high,k,h0,t1,t2,theta,a,b)$value
    Prop.death.treatment.split<-1-exp(-cum.hazard.treatment.split/split.interval.num)
    for(j in 1:split.interval.num)
    {
      if(t.high>tau-A)
      {
        Prob.censor=1/round(num+1-t.low*split.num-j)
      }
      else
      {
        Prob.censor=0
      }
      treatment.resd<-1-loss.treatment.split-Prop.death.treatment.split
      control.resd<-1-loss.control.split-Prop.death.control.split
      trans.matrix.uncensored.treatment<-rbind(c(1,0,0),c(0,1,0),c(loss.treatment.split,Prop.death.treatment.split,treatment.resd))
      trans.matrix.uncensored.control<-rbind(c(1,0,0),c(0,1,0),c(loss.control.split,Prop.death.control.split,control.resd))
      trans.matrix.admin.censor<-rbind(c(1,0,0),c(0,1,0),c(Prob.censor,0,1-Prob.censor))
      trans.matrix.treatment=trans.matrix.uncensored.treatment%*%trans.matrix.admin.censor
      trans.matrix.control=trans.matrix.uncensored.control%*%trans.matrix.admin.censor
      Trans.matrix.treatment=new("markovchain",states=c("Loss","Event","Surv.Trt"),transitionMatrix=trans.matrix.treatment,name="Survival")
      Trans.matrix.control=new("markovchain",states=c("Loss","Event","Surv.cont"),transitionMatrix=trans.matrix.control,name="Survival")
      i=(l-1)*split.interval.num+j
      Trans.matrix.list.control[[i]]=Trans.matrix.control
      Trans.matrix.list.treatment[[i]]=Trans.matrix.treatment
      Time=t.low+j*unit.interval.period
      Control.vector<-Control.matrix[i,2:4]
      Treatment.vector<-Treatment.matrix[i,2:4]
      Control.matrix[i+1,]<-c(Time,Control.vector*Trans.matrix.control)
      Treatment.matrix[i+1,]<-c(Time,Treatment.vector*Trans.matrix.treatment)
    }
  }
  Survival.markovchainList.control<- new("markovchainList", markovchains =Trans.matrix.list.control,name = "Survival Transition list")
  Survival.markovchainList.treatment<- new("markovchainList", markovchains =Trans.matrix.list.treatment,name = "Survival Transition list")
  Survival.markovchainList<-list(Survival.markovchainList.control=Survival.markovchainList.control,Survival.markovchainList.treatment=Survival.markovchainList.treatment)
  Control.matrix=data.frame(Control.matrix)
  Treatment.matrix=data.frame(Treatment.matrix)
  return(list(Control.matrix=Control.matrix,Treatment.matrix=Treatment.matrix,Survival.markovchainList=Survival.markovchainList))
}

marcov.survival.stat.distribution.chain.simple<-function(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)
  ## the function used to calculate the distribution parameters for the test statistics
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## split.interval.num: the splitted number of the smallest interval for an Markov chain
  ## interval.period: the duration of the smallest interval
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time(used for the calculation of Phi and the calculation of MERT)
  ## t1.true,t2.true:the minimum and maximum delay time(used in l(t) for the generation of the Markov chain)
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## Allc: the allocation ratio of the sample sizes of the treatment and control groups
  ## test.type:1:log-rank test;2:Maximin efficiency robustness test
  ## a,b:the shape parameter of the distribution of the delay time
  ## Phi.type:1:Phi obtained by numeric integration;2:estimated Phi
{
  list.survival.states.chain.simple<-marcov.survival.states.chain.simple(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1.true,t2.true,theta,a,b)
  ##obtain the proportions of the states at different time points by the Markov chain
  Pc<-1/(1+Allc)
  Pt<-Allc/(1+Allc)
  Control.matrix<-list.survival.states.chain.simple$Control.matrix
  Control.matrix[,2:4]<-Pc*Control.matrix[,2:4]
  Treatment.matrix<-list.survival.states.chain.simple$Treatment.matrix
  Treatment.matrix[,2:4]<-Pt*Treatment.matrix[,2:4]
  num=tau*split.interval.num/interval.period## the total number of the splited intervals
  ##Calc.matrix: the matrix is used to calculate the distribution parameters of the test statistics; it is made up of the following 16 elements
  Calc.matrix<-data.frame(matrix(0,nrow=num,ncol=21))
  colnames(Calc.matrix)<-c("d.cont","d.trt","d.all","Surv.cont.init","Surv.trt.init","Hazard.cont","Hazard.trt","Hazard.ratio","Surv.cont.end","Surv.trt.end","Cont.tot.exit","Trt.tot.exit","Ratio.Surv","Phi","Survival.time.init","Survival.time.end","Weights","d.porp","gamma","eta","omega")
  #d.cont: the death rate of the control group in a tiny interval
  #d.trt: the death rate of the treatment group in a tiny interval
  #d.all: the death rate in a tiny interval
  #Surv.cont: the proportion of the alive subjects in the control group
  #Surv.trt: the proportion of the alive subjects in the treatment group
  #Hazard.cont: the hazard rate of the control group
  #Hazard.trt: the hazard rate of the treatment group
  #Hazard.ratio: the ratio of the hazard rate of  the control group to that of the treatment group, also denoted by theta
  #Ratio.Surv: the ratio of the proportion of the alive subjects in the control group to that in the treatment group, also denoted by phi
  #Phi:a parameter which is used to calculate the weights of the Maximin Efficiency Robust Test
  #Survival.time: the survival time corresponding to the observation in the matrix
  #Weights:weights function
  #d.porp: the death proportion in the interval
  #gamma:phi*theta/(1+phi*theta)-phi/(1+phi)
  #eta:phi/(1+phi)^2
  #omega:phi*theta/(1+phi*theta)^2
  unit.interval<-interval.period/split.interval.num
  Death.cum=0
  if(test.type==1)
  {
    t1=0
    t2=0
  }
  if(t1==0)
  {
    Phi.t1=0
  }
  if(t2==0)
  {
    Phi.t2=0
  }
  if(Phi.type==1)
  {
    for(i in 1:num)
    {
      #Calc.matrix[i,1]=Control.matrix[i+1,3]-Control.matrix[i,3]
      Calc.matrix$d.cont[i]=Control.matrix$Event[i+1]-Control.matrix$Event[i]
      #Calc.matrix[i,2]=Treatment.matrix[i+1,3]-Treatment.matrix[i,3]
      Calc.matrix$d.trt[i]=Treatment.matrix$Event[i+1]-Treatment.matrix$Event[i]
      #Calc.matrix[i,3]=Calc.matrix[i,1]+Calc.matrix[i,2]
      Calc.matrix$d.all[i]=Calc.matrix$d.cont[i]+Calc.matrix$d.trt[i]
      Death.cum=Death.cum+Calc.matrix$d.all[i]
      #Calc.matrix[i,4]=Control.matrix[i,5]
      Calc.matrix$Surv.cont.init[i]=Control.matrix$Surv.cont[i]
      Calc.matrix$Surv.cont.end[i]=Control.matrix$Surv.cont[i+1]
      #Calc.matrix[i,5]=Treatment.matrix[i,5]
      Calc.matrix$Surv.trt.init[i]=Treatment.matrix$Surv.Trt[i]
      Calc.matrix$Surv.trt.end[i]=Treatment.matrix$Surv.Trt[i+1]
      #Calc.matrix[i,6]=-log(1-Calc.matrix[i,1]/Calc.matrix[i,4])/unit.interval
      Calc.matrix$Hazard.cont[i]=-log(1-Calc.matrix$d.cont[i]/Calc.matrix$Surv.cont.init[i])/unit.interval
      #Calc.matrix[i,7]=-log(1-Calc.matrix[i,2]/Calc.matrix[i,5])/unit.interval
      Calc.matrix$Hazard.trt[i]=-log(1-Calc.matrix$d.trt[i]/Calc.matrix$Surv.trt.init[i])/unit.interval
      #Calc.matrix[i,8]=Calc.matrix[i,6]/Calc.matrix[i,7]
      Calc.matrix$Hazard.ratio[i]=Calc.matrix$Hazard.cont[i]/Calc.matrix$Hazard.trt[i]
      ##Calc.matrix[i,9]=Calc.matrix[i,4]/Calc.matrix[i,5]
      Calc.matrix$Ratio.Surv[i]=Calc.matrix$Surv.cont.init[i]/Calc.matrix$Surv.trt.init[i]
      Calc.matrix$Trt.tot.exit[i]=Treatment.matrix$Loss[i+1]
      Calc.matrix$Cont.tot.exit[i]=Control.matrix$Loss[i+1]
      if(i==1)
      {
        integrand.low<-hazard.function.control(0,k,h0)*Pc*Pt
        integrand.high.numerator<-Calc.matrix$Hazard.cont[i]*Calc.matrix$Surv.cont.end[i]*(Pt-Calc.matrix$Trt.tot.exit[i])
        integrand.high.denominator<-(Pt-Calc.matrix$Trt.tot.exit[i])+(Pc-Calc.matrix$Cont.tot.exit[i])
        integrand.high<-integrand.high.numerator/integrand.high.denominator
        Calc.matrix$Phi[i]=(integrand.low+integrand.high)/2*unit.interval
      }
      else
      {
        integrand.low.numerator<-Calc.matrix$Hazard.cont[i-1]*Calc.matrix$Surv.cont.end[i-1]*(Pt-Calc.matrix$Trt.tot.exit[i-1])
        integrand.low.denominator<-(Pt-Calc.matrix$Trt.tot.exit[i-1])+(Pc-Calc.matrix$Cont.tot.exit[i-1])
        integrand.low<-integrand.low.numerator/integrand.low.denominator
        integrand.high.numerator<-Calc.matrix$Hazard.cont[i]*Calc.matrix$Surv.cont.end[i]*(Pt-Calc.matrix$Trt.tot.exit[i])
        integrand.high.denominator<-(Pt-Calc.matrix$Trt.tot.exit[i])+(Pc-Calc.matrix$Cont.tot.exit[i])
        integrand.high<-integrand.high.numerator/integrand.high.denominator
        Calc.matrix$Phi[i]=Calc.matrix$Phi[i-1]+(integrand.low+integrand.high)/2*unit.interval
      }
      Survival.time.init=Control.matrix[i,1]
      Survival.time.end=Control.matrix[i+1,1]
      Calc.matrix$Survival.time.init[i]=Survival.time.init
      Calc.matrix$Survival.time.end[i]=Survival.time.end
      ##remove the influences of machine errors, which results from the summation of the survival time
      if(round(Survival.time.end,8)==t1)
      {
        Phi.t1=Calc.matrix$Phi[i]
      }
      if(round(Survival.time.end,8)==t2)
      {
        Phi.t2=Calc.matrix$Phi[i]
      }
      #if(round(Survival.time.end,8)==tau)
      #{
      #  Phi.tau=Calc.matrix$Phi[i]
      #}
      if(i==num-1)
      {
        Phi.tau=Calc.matrix$Phi[i]
      }
    }
  }
  else
  {
    for(i in 1:num)
    {
      #Calc.matrix[i,1]=Control.matrix[i+1,3]-Control.matrix[i,3]
      Calc.matrix$d.cont[i]=Control.matrix$Event[i+1]-Control.matrix$Event[i]
      #Calc.matrix[i,2]=Treatment.matrix[i+1,3]-Treatment.matrix[i,3]
      Calc.matrix$d.trt[i]=Treatment.matrix$Event[i+1]-Treatment.matrix$Event[i]
      #Calc.matrix[i,3]=Calc.matrix[i,1]+Calc.matrix[i,2]
      Calc.matrix$d.all[i]=Calc.matrix$d.cont[i]+Calc.matrix$d.trt[i]
      Death.cum=Death.cum+Calc.matrix$d.all[i]
      #Calc.matrix[i,4]=Control.matrix[i,5]+Control.matrix[i,6]
      Calc.matrix$Surv.cont.init[i]=Control.matrix$Surv.cont[i]
      Calc.matrix$Surv.cont.end[i]=Control.matrix$Surv.cont[i+1]
      #Calc.matrix[i,5]=Treatment.matrix[i,5]+Treatment.matrix[i,6]
      Calc.matrix$Surv.trt.init[i]=Treatment.matrix$Surv.Trt[i]
      Calc.matrix$Surv.trt.end[i]=Treatment.matrix$Surv.Trt[i+1]
      #Calc.matrix[i,6]=-log(1-Calc.matrix[i,1]/Calc.matrix[i,4])/unit.interval
      Calc.matrix$Hazard.cont[i]=-log(1-Calc.matrix$d.cont[i]/Calc.matrix$Surv.cont.init[i])/unit.interval
      #Calc.matrix[i,7]=-log(1-Calc.matrix[i,2]/Calc.matrix[i,5])/unit.interval
      Calc.matrix$Hazard.trt[i]=-log(1-Calc.matrix$d.trt[i]/Calc.matrix$Surv.trt.init[i])/unit.interval
      #Calc.matrix[i,8]=Calc.matrix[i,6]/Calc.matrix[i,7]
      Calc.matrix$Hazard.ratio[i]=Calc.matrix$Hazard.cont[i]/Calc.matrix$Hazard.trt[i]
      ##Calc.matrix[i,9]=Calc.matrix[i,4]/Calc.matrix[i,5]
      Calc.matrix$Ratio.Surv[i]=Calc.matrix$Surv.cont.init[i]/Calc.matrix$Surv.trt.init[i]
      Calc.matrix$Trt.tot.exit[i]=Treatment.matrix$Loss[i+1]
      Calc.matrix$Cont.tot.exit[i]=Control.matrix$Loss[i+1]
      if(i==1)
      {
        #Calc.matrix[i,10]=Calc.matrix[i,3]*Calc.matrix[i,4]*Calc.matrix[i,5]/(Calc.matrix[i,4]+Calc.matrix[i,5])^2
        Calc.matrix$Phi[i]=Calc.matrix$d.all[i]*Calc.matrix$Surv.cont.init[i]*Calc.matrix$Surv.trt.init[i]/(Calc.matrix$Surv.cont.init[i]+Calc.matrix$Surv.trt.init[i])^2
      }
      else
      {
        #Calc.matrix[i,10]=Calc.matrix[i-1,10]+Calc.matrix[i,3]*Calc.matrix[i,4]*Calc.matrix[i,5]/(Calc.matrix[i,4]+Calc.matrix[i,5])^2
        Calc.matrix$Phi[i]=Calc.matrix$Phi[i-1]+Calc.matrix$d.all[i]*Calc.matrix$Surv.cont.init[i]*Calc.matrix$Surv.trt.init[i]/(Calc.matrix$Surv.cont.init[i]+Calc.matrix$Surv.trt.init[i])^2
      }
      Survival.time.init=Control.matrix[i,1]
      Survival.time.end=Control.matrix[i+1,1]
      Calc.matrix$Survival.time.init[i]=Survival.time.init
      Calc.matrix$Survival.time.end[i]=Survival.time.end
      ##remove the influences of machine errors, which results from the summation of the survival time
      if(round(Survival.time.end,8)==t1)
      {
        Phi.t1=Calc.matrix$Phi[i]
      }
      if(round(Survival.time.end,8)==t2)
      {
        Phi.t2=Calc.matrix$Phi[i]
      }
      if(round(Survival.time.end,8)==tau)
      {
        Phi.tau=Calc.matrix$Phi[i]
      }
    }
  }
  for(i in 1:num)
  {
    Survival.time.end=Calc.matrix$Survival.time.end[i]
    ##Caution:<=t1 or <t1,please determine by debug
    if(round(Survival.time.end,8)<t1)
    {
      Calc.matrix$Weights[i]=0
    }
    else if(round(Survival.time.end,8)>=t1 & round(Survival.time.end,8)<=t2)
    {
      ##Caution:>=t1 or >t1,please determine by debug
      Calc.matrix$Weights[i]=((Phi.tau-Calc.matrix$Phi[i])/(Phi.tau-Phi.t1))^(-0.5)
      #Note:tau is commonly assumed to be greated than t2; therefore, don't worry that the above value is Inf
    }
    else
    {
      Calc.matrix$Weights[i]=2*((Phi.tau-Phi.t2)/(Phi.tau-Phi.t1))^(-0.5)
    }
    Calc.matrix$d.porp[i]=Calc.matrix$d.all[i]/Death.cum
  }
  ##The parameters for the number of deaths
  Mean.unstd.d.unit=0
  Var.H0.unstd.d.unit=0
  Var.H1.unstd.d.unit=0
  if(test.type==1)
  {
    Calc.matrix$Weights=Calc.matrix$Weights/2
    ##log-rank test
    for(i in 1:num)
    {
      phi=Calc.matrix$Ratio.Surv[i]
      theta=Calc.matrix$Hazard.ratio[i]
      gamma=phi*theta/(1+phi*theta)-phi/(1+phi)
      #Calc.matrix[i,14]=gamma
      Calc.matrix$gamma[i]=gamma
      eta=phi/(1+phi)^2
      #Calc.matrix[i,15]=eta
      Calc.matrix$eta[i]=eta
      omega=phi*theta/(1+phi*theta)^2
      #Calc.matrix[i,16]=omega
      Calc.matrix$omega[i]=omega
      #rho=Calc.matrix[i,13]
      rho=Calc.matrix$d.porp[i]
      Mean.unstd.d.unit=Mean.unstd.d.unit+gamma*rho
      Var.H0.unstd.d.unit=Var.H0.unstd.d.unit+eta*rho
      Var.H1.unstd.d.unit=Var.H1.unstd.d.unit+omega*rho
    }
  }
  else
  {
    ## Maximin efficiency robust test
    for(i in 1:num)
    {
      weight=Calc.matrix$Weights[i]
      phi=Calc.matrix$Ratio.Surv[i]
      theta=Calc.matrix$Hazard.ratio[i]
      gamma=phi*theta/(1+phi*theta)-phi/(1+phi)
      #Calc.matrix[i,14]=gamma
      Calc.matrix$gamma[i]=gamma
      eta=phi/(1+phi)^2
      #Calc.matrix[i,15]=eta
      Calc.matrix$eta[i]=eta
      omega=phi*theta/(1+phi*theta)^2
      #Calc.matrix[i,16]=omega
      Calc.matrix$omega[i]=omega
      #rho=Calc.matrix[i,13]
      rho=Calc.matrix$d.porp[i]
      Mean.unstd.d.unit=Mean.unstd.d.unit+gamma*rho*weight
      Var.H0.unstd.d.unit=Var.H0.unstd.d.unit+rho*eta*weight^2
      Var.H1.unstd.d.unit=Var.H1.unstd.d.unit+rho*omega*weight^2
    }
  }
  Mean.unstd.d.unit=as.numeric(Mean.unstd.d.unit)
  Var.H0.unstd.d.unit=as.numeric(Var.H0.unstd.d.unit)
  Var.H1.unstd.d.unit=as.numeric(Var.H1.unstd.d.unit)
  return(list(Mean.unstd.d.unit=Mean.unstd.d.unit,Var.H0.unstd.d.unit=Var.H0.unstd.d.unit,Var.H1.unstd.d.unit=Var.H1.unstd.d.unit,Calc.matrix=Calc.matrix,Death.cum=as.numeric(Death.cum),list.survival.states.chain.simple=list.survival.states.chain.simple))
}

death.estimation.fixed.sample.design<-function(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type,alpha,power,sides)
  ## the function used to calculate the distribution parameters for the test statistics
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## split.interval.num: the splitted number of the smallest interval for an Markov chain
  ## interval.period: the duration of the smallest interval
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## noncomp.control: the proportion of the non-compliance subjects in the control group in an unit period
  ## noncomp.treatment: the proportion of the non-compliance subjects in the treatment group in an unit period
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2: the minimum and maximum delay time(which is used in the MERT)
  ## t1.true,t2.true:the minimum and maximum delay time(which is used in the delay function l(t) and the generation of the Markov chain)
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
  ## Allc: the allocation ratio of the sample sizes of the treatment and control groups
  ## test.type:1:log-rank test;2:Maximin efficiency robustness test
  ## alpha: the significance level
  ## power: the nominal power
  ## sides: 1:right-sided test;2:two-sided test
  ## Phi.type:1:Phi obtained by numeric integration;2:estimated Phi
{
  list.distribution.chain.simple<-marcov.survival.stat.distribution.chain.simple(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)
  Mean.unstd.d.unit<-list.distribution.chain.simple$Mean.unstd.d.unit
  Var.H0.unstd.d.unit<-list.distribution.chain.simple$Var.H0.unstd.d.unit
  Var.H1.unstd.d.unit<-list.distribution.chain.simple$Var.H1.unstd.d.unit
  if(sides==2)
  {
    d=((qnorm(1-alpha/2)*sqrt(Var.H0.unstd.d.unit)+qnorm(power)*sqrt(Var.H1.unstd.d.unit))/Mean.unstd.d.unit)^2
  }
  else
  {
    d=((qnorm(1-alpha)*sqrt(Var.H0.unstd.d.unit)+qnorm(power)*sqrt(Var.H1.unstd.d.unit))/Mean.unstd.d.unit)^2
  }
  death.proportion<-list.distribution.chain.simple$Death.cum
  n=round(d/death.proportion)
  nc=ceiling(n/(1+Allc))
  nt=n-nc
  d=n*death.proportion
  Mean=sqrt(d)*Mean.unstd.d.unit/sqrt(Var.H0.unstd.d.unit)
  Var=Var.H1.unstd.d.unit/Var.H0.unstd.d.unit
  Asymptotic.distribution<-list(Mean=Mean,Var=Var)
  return(list(d=d,n=n,nc=nc,nt=nt,list.distribution.chain.simple=list.distribution.chain.simple,Asymptotic.distribution=Asymptotic.distribution))
}

death.estimation.fixed.sample.design.explore<-function(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type,alpha,power,sides)
  ## the function used to calculate the distribution parameters for the test statistics(the variance is assumed to be the squart root of the original one)
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## split.interval.num: the splitted number of the smallest interval for an Markov chain
  ## interval.period: the duration of the smallest interval
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## noncomp.control: the proportion of the non-compliance subjects in the control group in an unit period
  ## noncomp.treatment: the proportion of the non-compliance subjects in the treatment group in an unit period
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2: the minimum and maximum delay time(which is used in the MERT)
  ## t1.true,t2.true:the minimum and maximum delay time(which is used in the delay function l(t) and the generation of the Markov chain)
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
  ## Allc: the allocation ratio of the sample sizes of the treatment and control groups
  ## test.type:1:log-rank test;2:Maximin efficiency robustness test
  ## alpha: the significance level
  ## power: the nominal power
  ## sides: 1:right-sided test;2:two-sided test
  ## Phi.type:1:Phi obtained by numeric integration;2:estimated Phi
{
  list.distribution.chain.simple<-marcov.survival.stat.distribution.chain.simple(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)
  Mean.unstd.d.unit<-list.distribution.chain.simple$Mean.unstd.d.unit
  Var.H0.unstd.d.unit<-list.distribution.chain.simple$Var.H0.unstd.d.unit
  Var.H1.unstd.d.unit<-list.distribution.chain.simple$Var.H1.unstd.d.unit
  if(sides==2)
  {
    d=((qnorm(1-alpha/2)*sqrt(Var.H0.unstd.d.unit)+qnorm(power)*(Var.H0.unstd.d.unit*Var.H1.unstd.d.unit)^0.25)/Mean.unstd.d.unit)^2
  }
  else
  {
    d=((qnorm(1-alpha)*sqrt(Var.H0.unstd.d.unit)+qnorm(power)*(Var.H0.unstd.d.unit*Var.H1.unstd.d.unit)^0.25)/Mean.unstd.d.unit)^2
  }
  death.proportion<-list.distribution.chain.simple$Death.cum
  n=round(d/death.proportion)
  nc=ceiling(n/(1+Allc))
  nt=n-nc
  d=n*death.proportion
  Mean=sqrt(d)*Mean.unstd.d.unit/sqrt(Var.H0.unstd.d.unit)
  Var=sqrt(Var.H1.unstd.d.unit/Var.H0.unstd.d.unit)
  Asymptotic.distribution<-list(Mean=Mean,Var=Var)
  return(list(d=d,n=n,nc=nc,nt=nt,list.distribution.chain.simple=list.distribution.chain.simple,Asymptotic.distribution=Asymptotic.distribution))
}

death.estimation.fixed.sample.design.old<-function(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type,alpha,power,sides)
  ## the function used to calculate the distribution parameters for the test statistics
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## split.interval.num: the splitted number of the smallest interval for an Markov chain
  ## interval.period: the duration of the smallest interval
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## noncomp.control: the proportion of the non-compliance subjects in the control group in an unit period
  ## noncomp.treatment: the proportion of the non-compliance subjects in the treatment group in an unit period
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2: the minimum and maximum delay time(which is used in the MERT)
## t1.true,t2.true:the minimum and maximum delay time(which is used in the delay function l(t) and the generation of the Markov chain)
## theta: the hazard ratio after complete onset of the treatment effect
## a,b:the shape parameter of the distribution of the delay time
## Allc: the allocation ratio of the sample sizes of the treatment and control groups
## test.type:1:log-rank test;2:Maximin efficiency robustness test
## alpha: the significance level
## power: the nominal power
## sides: 1:right-sided test;2:two-sided test
## Phi.type:1:Phi obtained by numeric integration;2:estimated Phi
{
  list.distribution.chain.simple<-marcov.survival.stat.distribution.chain.simple(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)
  Mean.unstd.d.unit<-list.distribution.chain.simple$Mean.unstd.d.unit
  Var.H0.unstd.d.unit<-list.distribution.chain.simple$Var.H0.unstd.d.unit
  if(sides==2)
  {
    d=Var.H0.unstd.d.unit*((qnorm(1-alpha/2)+qnorm(power))/Mean.unstd.d.unit)^2
  }
  else
  {
    d=Var.H0.unstd.d.unit*((qnorm(1-alpha)+qnorm(power))/Mean.unstd.d.unit)^2
  }
  death.proportion<-list.distribution.chain.simple$Death.cum
  n=round(d/death.proportion)
  nc=ceiling(n/(1+Allc))
  nt=n-nc
  d=n*death.proportion
  Mean=sqrt(d)*Mean.unstd.d.unit/sqrt(Var.H0.unstd.d.unit)
  Var=1
  Asymptotic.distribution<-list(Mean=Mean,Var=Var)
  return(list(d=d,n=n,nc=nc,nt=nt,list.distribution.chain.simple=list.distribution.chain.simple,Asymptotic.distribution=Asymptotic.distribution))
}

