###----------------------------------------------------------------###
# Project Name : Count Models
# Code Name: R - Project2
###----------------------------------------------------------------###

###----------------------------------------------------------------###
### Contents
## 1   Poisson Model
## 1.1 Hurdle Poisson Model
## 1.2 Zero Inflated Poisson Model
## 2   Negative Binomial Model
## 2.1 Hurdle Negative Binomial Model
## 2.2 Zero Inflated Negative Binomial Model
## 3   Model Evaluation
###----------------------------------------------------------------###

##Installing and Loading Libraries
install.packages("pscl")
library(pscl)
install.packages("lmtest")
library(lmtest)
install.packages("bbmle")
library(bbmle)
###----------------------------------------------------------------###
#Import the Health dataset
Health <- read.csv("C:DUMMY/Health.csv")
View(Health)


###----------------------------------------------------------------###
#Preliminary data analysis
#Histogram of health visits
plot(table(Health$ofp),col="darkgrey")

#The above plot indicates that there is high variation in the hospital visits(Over dispersion problem)
#and also high number of the zeros in the data (Excess zeros probelm) 
#The project scope is limited to handling excess zero's problem using hurdle and Zero inflated models

###----------------------------------------------------------------###
###--------------------1. Poisson Models---------------------------###
###----------------------------------------------------------------###

# Poisson Regression using glm
fn_poisson = glm(Health$ofp~ Health$numchron+Health$married+Health$employed,data = Health, family = poisson)
#Excess zeros is not captured by this classical model
summary(fn_poisson)


# Poisson Regression using Maximaization of Loglikelihood function
LLpoisson = function(b0,b1,b2,b3)
{
  y= b0 + b1*(Health$numchron) + b2*(Health$employed) + b3*(Health$married)
  l= exp(y)
  
  LL=sum(dpois(Health$ofp,lambda = l,log = T))
  return(-1*LL)
}

# Run the MLE2 function by giving list of start values of the parameters 
# For maximum likelihood, this function returns the parameter estimates

# Refer document for explaination of  start values
res1 = mle2(minuslogl = LLpoisson,start = list(b0=log(mean(Health$ofp),base = exp(1)), b1=0, b2=0, b3=0))
summary(res1)

###----------------------------------------------------------------###
###----------------1.1 Hurdle Poisson Model------------------------###
###----------------------------------------------------------------###

# Handling Over dispersion and excess zeros problem in count data using hurdle models
fn_hurdle_pois<- hurdle(ofp ~ ., data = Health, dist = "poisson")
summary(fn_hurdle_pois)

# Hurdle model with poisson using Likelihood function
LLHurdle = function(a0,a1,a2,b0,b1,b2,b3)
{
  x = a0 + a1*Health$numchron + a2*Health$male
  p = exp(x)/(1+exp(x)) # transforming x to get the value in range 0 to 1
  
  y = b0 + b1*Health$numchron + b2*Health$employed + b3*Health$married
  lam = exp(y) # transforming y to get non negative values
  
  L = ifelse(Health$ofp == 0,p,
      ((1-p)/(1-dpois(0,lambda = lam)))*dpois(Health$ofp,lambda = lam))
  LLsum = sum(log(L))
  return(-1*LLsum)
}

# Run the MLE2 function by giving list of start values of the parameters 
# For maximum likelihood, this function returns the parameter estimates

# Refer document for explaination of  start values
res4 = mle2(minuslogl = LLHurdle,start = list(a0=0,a1=0,a2=0,b0=log(mean(Health$ofp),base = exp(1)),b1=0,b2=0,b3=0))
summary(res4)



###----------------------------------------------------------------###
###------------1.2 Zero Inflated Poisson Model---------------------###
###----------------------------------------------------------------###

# Handling excess zeros problem in count data using Zero inflated models
fn_zi_pois <- zeroinfl(ofp ~ + Health$numchron+Health$employed+Health$married|Health$numchron+Health$male, data = Health, dist = "poisson")
summary(fn_zi_pois)


# Zero Inflated with Poisson dist.
LLZeroPois = function(a0,a1,a2,b0,b1,b2,b3)
{
  x = a0 + a1*Health$numchron + a2*Health$male
  p = exp(x)/(1+exp(x)) # transforming x to get the value in range 0 to 1
  
  y = b0 + b1*Health$numchron + b2*Health$employed + b3*Health$married
  lam = exp(y) # transforming y to get non negative values
  
  L = ifelse(Health$ofp == 0, p+((1-p)*(dpois(0,lambda = lam))), (1-p)*(dpois(Health$ofp,lambda = lam)))
  LLsum = sum(log(L))
  return(-1*LLsum)
}

# Run the MLE2 function by giving list of start values of the parameters 
# For maximum likelihood, this function returns the parameter estimates

# Refer document for explaination of  start values
res2 = mle2(minuslogl = LLZeroPois,start = list(a0=0,a1=0,a2=0,b0=log(mean(Health$ofp),base = exp(1)),b1=0,b2=0,b3=0))
summary(res2)



###----------------------------------------------------------------###
###--------------2 Negative Binomial Models------------------------###
###----------------------------------------------------------------###

#Negative binomial regression using glm
fn_nbin = glm.nb(ofp ~ ., data = Health)
summary(fn_nbin)


# Negative Bunimial using Maximaization of Loglikelihood function
LLNB = function(b0,b1,b2,b3)
{ 
  y= b0 + b1*(Health$numchron) + b2*(Health$employed) + b3*(Health$married)
  l= exp(y)
  
  prob=dnbinom(Health$ofp, size=1000 , mu= l , log = TRUE)
  LLsum = sum(prob)
  return(-1*LLsum) 
}

# Run the MLE2 function by giving list of start values of the parameters 
# For maximum likelihood, this function returns the parameter estimates

# Refer document for explaination of  start values

res2 = mle2(minuslogl = LLNB, start = list( b0=log(mean(Health$ofp),base = exp(1)), 
                                           b1=0, b2=0, b3=0))
summary(res2)

###----------------------------------------------------------------###
###----------2.1 Hurdle Negative Binomial Model--------------------###
###----------------------------------------------------------------###

# Handling Over dispersion and excess zeros problem in count data using hurdle models
fn_hurdle_nb <- hurdle(ofp ~ ., data = Health, dist = "negbin")
summary(fn_hurdle_nb)


# Hurdle model with Negative binomial 
LLNbHurdle = function(a0,a1,a2,b0,b1,b2,b3)
{
  x = a0 + a1*Health$numchron + a2*Health$male
  p = exp(x)/(1+exp(x)) # transforming x to get the value in range 0 to 1
  
  y = b0 + b1*Health$numchron + b2*Health$employed + b3*Health$married
  lam = exp(y) # transforming y to get non negative values
  
  L = sum(ifelse(Health$ofp==0,p,((1-p)/(1-dnbinom(0,size =1000,mu = lam)))*dnbinom(Health$ofp,size =1000,mu = lam)))
  LLsum = sum(log(L))
  return(-1*LLsum)
}

# Run the MLE2 function by giving list of start values of the parameters 
# For maximum likelihood, this function returns the parameter estimates

# Refer document for explaination of  start values
res1 = mle2(minuslogl = LLNbHurdle,start = list(a0=0,a1=0,a2=0,b0=log(mean(Health$ofp),base = exp(1)),b1=0,b2=0,b3=0))
summary(res1)

###----------------------------------------------------------------###
###----------2.2 Zero Inflated Negative Binomial Model-------------###
###----------------------------------------------------------------###
   
# Using zeroinfl function to develop zero inflated negative binomial model
fn_zinb <- zeroinfl(ofp ~ + Health$numchron+Health$employed+Health$married|Health$numchron+Health$male, data = Health, dist = "negbin")
summary(fn_zinb)


# Zero Inflated with Negative Binomial by maximizing loglikelihood function
LLZeroNb = function(a0,a1,a2,b0,b1,b2,b3)
{
  x = a0 + a1*Health$numchron + a2*Health$male
  p = exp(x)/(1+exp(x)) # transforming x to get the value in range 0 to 1
  
  y = b0 + b1*Health$numchron + b2*Health$employed + b3*Health$married
  lam = exp(y) # transforming y to get non negative values
  
  L = ifelse(Health$ofp == 0, p+((1-p)*(dnbinom(0,size =1000,mu = lam))), (1-p)*(dnbinom(Health$ofp,size =1000,mu = lam)))
  LLsum = sum(log(L))
  return(-1*LLsum)
}

# Run the MLE2 function by giving list of start values of the parameters 
# For maximum likelihood, this function returns the parameter estimates

# Refer document for explaination of  start values
res3 = mle2(minuslogl = LLZeroNb,start = list(a0=0,a1=0,a2=0,b0=log(mean(Health$ofp),base = exp(1)),b1=0,b2=0,b3=0))
summary(res3)


###----------------------------------------------------------------###
###----------------------3 Model Evaluation------------------------###
###----------------------------------------------------------------###

AIC(fn_poisson, fn_hurdle_pois, fn_zi_pois, fn_nbin, fn_hurdle_nb, fn_zinb)

fk=logLik(fn_poisson)
logLik(fn_hurdle_pois)
logLik(fn_zi_pois)
logLik(fn_nbin)
logLik(fn_hurdle_nb)
logLik(fn_zinb)


#Creating a new dataframe named "Results" for comparing model results

#Append AIC values of models

results=AIC(fn_poisson, fn_hurdle_pois, fn_zi_pois, fn_nbin, fn_hurdle_nb, fn_zinb)

rownames(results) = NULL

#Append loglikelihood values of all models
results$LogLik=rbind(logLik(fn_poisson),
                     logLik(fn_hurdle_pois),
                     logLik(fn_zi_pois),
                     logLik(fn_nbin),
                     logLik(fn_hurdle_nb),
                     logLik(fn_zinb))

#Specifying model names to the results dataframe 
results$models = rbind("fn_poisson", "fn_hurdle_pois", "fn_zi_pois", "fn_nbin", "fn_hurdle_nb", "fn_zinb")

results=results[,c(4,1,2,3)]

#Viewing the model results
View(results)
