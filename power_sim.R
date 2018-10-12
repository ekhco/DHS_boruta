
# ------------------------------------------------------------------------------
# power simulation for logistic regression

lin_int = function(n, se) {  # for n participatns
  m = 3*n  # because 3 observations per person
  CBT = rep(0:1,each=round(m/2,0))  # the treatment, half and half       THE FACTOR
  x = sample(c(0,1,30),length(CBT), replace = TRUE)           # normally distributed random x   THE CONTINUOUS (time?)
  table(CBT,x)
  e = rnorm(length(CBT), 3.4, se*sqrt(n))
  lp = e + 0.028*CBT - 0.00197*x - 6*0.00145*CBT*x        # this is the interaction!
  y = lp
  data <- data.frame(y,x,CBT=factor(CBT))
  head(data)


  # have to do link things for log to transform into a 1/0
  # link_lp = exp(lp)/(1 + exp(lp))
  # y = (runif(n) < lp)

  # this is the linear model w interaction
  lin.int = glm(y ~ CBT*x, data = data)
  summary(lin.int)
  reject = ifelse( coef(summary(lin.int))[4,4] < .05, 1, 0)
      # The coef() function above gets the parameter estimates; the [4,4]
      # element is the p-value for the interaction.
  return(reject)
}
# Running the function many times is also trivial, using the replicate() function.

pow1 = replicate(199, lin_int(n=142, se = 0.0315)) # 0.0063))
#Get the estimated power and confidence limits, we use the binom.test() function.
binom.test(sum(pow1), 199)


# histogram of STAI
n = 142
sd = 0.0315*sqrt(n)
e = rnorm(3*n, 3.4, sd)
hist(exp(e))













# power simulation for logistic regression

logist_inter = function(n) {
  c = rep(0:1,each=n/2)  # sample size is 100 THE FACTOR
  x = rnorm(n)           # normally distributed random x
  lp = -3 + 2*c*x        # this is the interaction!
  # have to do link things for log to transform into a 1/0
  link_lp = exp(lp)/(1 + exp(lp))
  y = (runif(n) < link_lp)

  log.int = glm(y~as.factor(c)*x, family=binomial)
  reject = ifelse( coef(summary(log.int))[4,4] < .05, 1, 0)
      # The coef() function above gets the parameter estimates; the [4,4]
      # element is the p-value for the interaction.
  return(reject)
}
# Running the function many times is also trivial, using the replicate() function.

pow1 = replicate(199, logist_inter(200))
# The result is an array of 1s and 0s. To get the estimated power and confidence limits, we use the binom.test() function.

binom.test(sum(pow1), 199)
