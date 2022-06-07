##
## Estimating JS and KL divergences between sIT models
##

## Method of moments estimate of the beta distribution parameters from the data
alpha.beta.est <- function(pp)
{
    m <- mean(pp)
    v <- var(pp)

    m1m <- m*(1-m)

    if ( v > m1m )
        stop(paste0("v > m*(1-m): v=",v," m*(1-m)=",m1m))

    z <- m1m/v - 1
    a <- m*z
    b <- (1-m)*z

    list(a=a, b=b, v=v, m1m=m1m)
}

source("~/devel/speciateIT/R/rLpps_vs_observed_llps/plot_dlpp_lib.R")

pics.dir <- "~/devel/speciateIT/R/JS_KL_bw_sIT_models/pics/"

file2 <- "~/devel/speciateIT/data/vaginal_Mar2021_V3V4/V3V4_vag_mcDir/rLpps/BVAB1.csv"
rlpps <- read.csv(file2, header=FALSE, row.names=1)

(ref.tx <- rownames(rlpps)[1])

rr <- plot.drlpp(rlpps, ref.tx)

pp <- 10^as.numeric(rlpps[1,])
myHist(pp, freq=F, xlim=c(0,1))


r <- alpha.beta.est(pp)
x <- seq(0,1, length=500)
db.mm <- dbeta(x, shape1=r$a, shape2=r$b)
myHist(pp, freq=F)
lines(x, db.mm, col='red')

## Finding an optimal value of the shape par's using optim()

h <- myHist(pp, freq=F)
d <- density(pp, from=0, to=1)
lines(d$x, d$y)
lines(h$mids, h$density)

opt.fn <- function(x)
{
    db <- dbeta(h$mids, shape1=x[1], shape2=x[2])
    sum((db - h$density)^2)
}



## stan model

opt_chr1 <- "
data {
  int<lower=0> N;
  real x[N];
}
parameters {
  real<lower = 0> alpha0;
  real<lower = 0> beta0;
}
model {
  alpha0 ~ normal(0, 1);
  beta0 ~ normal(0, 10);
  //target += beta_lpdf(x | alpha0, beta0); // same as below
  x ~ beta(alpha0, beta0);
}
"

library(rstan)
## compile model (~ 10 to 15 seconds)
##opt_mod1 <- stan_model(model_code = opt_chr1)

## initialize parameter values (based on knowledge or fitdist results)
init_list <- list(alpha0 = 7, beta0 = 1)

# optimize data given model
opt1 <- optimizing(object = opt_mod1, as_vector = FALSE,
                  data = list(x = pp,
                              N = length(pp)),
                  hessian = TRUE,
                  draws = 2500)
opt1$par
x <- seq(0,1, length=500)
db.stan <- dbeta(x, shape1=opt1$par$alpha0, shape2=opt1$par$beta0)
myHist(pp, freq=F)
lines(x, db.stan, col='blue')
lines(x, db.mm, col='red')

db.c <- dbeta(x, shape1=15, shape2=1.1)
lines(x, db.c, col='green')

## =============================================
##
##   JS estimates using mm estimates of beta
##
## =============================================

pp.ref <- 10^as.numeric(rlpps[1,])
pp.sib <- 10^as.numeric(rlpps[2,])

dx <- 0.01
x <- seq(dx,1-dx, length=500)

r <- alpha.beta.est(pp.ref)
db.ref <- dbeta(x, shape1=r$a, shape2=r$b)

r <- alpha.beta.est(pp.sib)
db.sib <- dbeta(x, shape1=r$a, shape2=r$b)

db.m <- (db.ref + db.sib) / 2

kl.div <- function(d1, d2)
{
    idx <- d1>0 & d2>0
    d1 <- d1[idx]
    d2 <- d2[idx]

    mean(d1*log(d1/d2))
}

kl.div(db.ref, db.m)

0.5 * ( kl.div(db.ref, db.m) + kl.div(db.sib, db.m) ) # 0.6754181

kl.fn <- approxfun(x, y=db.ref*log(db.ref/db.m))
integrate(kl.fn, lower=dx, upper=(1-dx)) # 0.6283301
kl.div(db.ref, db.m)                     # 0.6455245

kl.fn <- approxfun(x, y=db.sib*log(db.sib/db.m))
integrate(kl.fn, lower=dx, upper=(1-dx)) # 0.6925799
kl.div(db.sib, db.m)                     # 0.7053117


## =============================================
##
##   JS estimates using stan estimates of beta
##
## =============================================

pp.ref <- 10^as.numeric(rlpps[1,])
pp.sib <- 10^as.numeric(rlpps[2,])

m.ref <- optimizing(object = opt_mod1, as_vector = FALSE,
                   data = list(x = pp.ref,
                               N = length(pp)),
                   hessian = TRUE,
                   draws = 2500)
m.ref$par
dx <- 0.01
x <- seq(dx,1-dx, length=500)
db.ref.stan <- dbeta(x, shape1=m.ref$par$alpha0, shape2=m.ref$par$beta0)


m.sib <- optimizing(object = opt_mod1, as_vector = FALSE,
                   data = list(x = pp.sib,
                               N = length(pp)),
                   hessian = TRUE,
                   draws = 2500)
m.sib$par
dx <- 0.01
x <- seq(dx,1-dx, length=500)
db.sib.stan <- dbeta(x, shape1=m.sib$par$alpha0, shape2=m.sib$par$beta0)

db.m.stan <- (db.ref.stan + db.sib.stan) / 2

0.5 * ( kl.div(db.ref.stan, db.m.stan) + kl.div(db.sib.stan, db.m.stan) ) # 0.6868376

kl.fn <- approxfun(x, y=db.ref.stan*log(db.ref.stan/db.m.stan))
integrate(kl.fn, lower=dx, upper=(1-dx)) # 0.6531718
kl.div(db.ref.stan, db.m.stan)           # 0.6692877

kl.fn <- approxfun(x, y=db.sib.stan*log(db.sib.stan/db.m.stan))
integrate(kl.fn, lower=dx, upper=(1-dx)) # 0.6916813
kl.div(db.sib.stan, db.m.stan)           # 0.7043876



##
## JS of ref and sib
##

#'
#'
js.fn <- function(i, rlpps, dx=0.01)
{
    pp.ref <- 10^as.numeric(rlpps[1,])
    pp.sib <- 10^as.numeric(rlpps[i,])

    x <- seq(dx,1-dx, length=500)

    r <- alpha.beta.est(pp.ref)
    db.ref <- dbeta(x, shape1=r$a, shape2=r$b)

    r <- alpha.beta.est(pp.sib)
    db.sib <- dbeta(x, shape1=r$a, shape2=r$b)

    db.m <- (db.ref + db.sib) / 2

    0.5 * ( kl.div(db.ref, db.m) + kl.div(db.sib, db.m) )
}

js.fn(4, rlpps)


##
## computing JS
##

dir <- "~/devel/speciateIT/data/vaginal_Mar2021_V3V4/V3V4_vag_mcDir/rLpps/"
files <- dir(dir, pattern =".csv")

js.est <- c()
counter <- 1
for ( file in files )
{
    cat("\r", counter)
    counter <- counter + 1
    file <- paste0(dir,file)
    rlpps <- read.csv(file, header=FALSE, row.names=1)
    nr <- nrow(rlpps)
    for ( i in 2:nr ) {
        js.est <- c(js.est, js.fn(i, rlpps))
    }
}

myHist(js.est)
