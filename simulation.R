############################################################################
#
# Simulate indirect and direct effects with right and interval
#   censoring
# created: May 03, 2022
#
#
############################################################################


# load libraries ----------------------------------------------------------

library(readxl)
library(tidyverse)
library(survival)
library(interval)

# source necessary functions --------------------------------------------------------

source("functions/rmst_ic.R")
source("functions/pseudo_rmst.R")

# set parameter values ----------------------------------------------------

nsample <- 100

tstar <- 8
lam <- 1.5
nu <- 0.8

# simulate covariates -----------------------------------------------------

Csim <- rbinom(nsample, 1, 0.5)
Asim <- rbinom(nsample, 1, 0.5)

eta <- -0.10 + 1.09 * Asim - 0.43 * Csim
psim <- 1/(1+exp(-eta))

Msim <- rbinom(nsample, 1, psim)

# simulate time to events -------------------------------------------------

X <- cbind(Asim, Msim, Csim)

Xbeta <- X %*% c(-0.2, -1, -0.3)

#sample N from a uniform(0, 1) distribution

U <- stats::runif(nsample)

#calculate simulated times

t <- c((-(log(U)) / (lam * exp(Xbeta)))^(1/nu))

# simulate right censoring ------------------------------------------------

#create censoring variable
ct <- runif(nsample, 0, 100)

t_rc <- ifelse(t < ct, t, ct)

event_rc <- t < ct

# simulate intervals  ---------------------------------------------------

start <- runif(nsample, 0, 0.5)

# matrix of study visit times ---------------------------------------------

tmat <- matrix(start, 
               ncol = 21, 
               nrow = nsample) +
  matrix((0:20) * 0.5, 
         byrow = T, 
         ncol = 21, 
         nrow = nsample)


# simulate missed visits 

miss <- cbind(1, matrix(rbinom(20*nsample, 1, 0.8), 
                        nrow = nsample, 
                        ncol = 20))

tmat_miss <- miss*tmat

# create intervals 

L <- tmat_miss *(tmat_miss < t)
l <- apply(L, 1, function(x)max(x, na.rm = T))
R <- tmat_miss * (tmat_miss >= t)
R[R == 0] <- NA          
r <- apply(R, 1,function(x)min(x, na.rm = T))


# simulated dataframe -----------------------------------------------------

#output as data.frame
tdata <- data.frame(t = t,
                    t_rc = t_rc,
                    event_rc = event_rc,
                    l, r, 
                    tmid = (l + r)/2,
                    event_mid = as.numeric((l + r)/2 < Inf),
                    X)

# right censoring ----------------------------------------------------------------

#RMST for A = 0 

sfit_strat_rc <- survfit(Surv(tdata$t_rc, 
                              event = tdata$event_rc) ~ tdata$Asim)

rmsts_rc <- summary(sfit_strat_rc, rmean = tstar)$table

theta0_rc <- rmsts_rc[1 ,"rmean"]
theta1_rc <- rmsts_rc[2 ,"rmean"]

total_rc <- theta1_rc - theta0_rc

#outcome model 

sfit_rc <- survfit(Surv(tdata$t_rc, 
                        event = tdata$event_rc) ~ 1)

tdata$psvals_rc <- pseudo(sfit_rc, times = 8, type = "rmst")

#fit outcome model on pseduo-observations

y_model_rc <- lm(psvals_rc ~ Msim + Csim,
                 subset = which(Asim == 0),
                 data = tdata)

#estimate E[Y^(0, I = 1)]

ey0i_rc <- predict(y_model_rc, 
                   newdata = tdata[Asim == 1,])


#indirect and direct effects

ie_rc <- mean(ey0i_rc) - theta0_rc
de_rc <- theta1_rc - mean(ey0i_rc)


# interval censoring ------------------------------------------------------

# RMST total effect 

theta0_ic <- rmst_ic(l, r,
                     subset = (Asim == 0),
                     tau = tstar)
theta1_ic <- rmst_ic(l, r,
                     subset = (Asim == 1),
                     tau = tstar)

total_ic <- theta1_ic - theta0_ic

# outcome model -----------------------------------------------------------

#calculate pseudo-observations
pseudo_ic <- pseudo_rmst(left = l[Asim == 0], right = r[Asim == 0], 
                         tau = tstar, 
                         censor = "interval")

Msim0 <- Msim[Asim == 0]
Csim0 <- Csim[Asim == 0]

#fit model on pseduo-observations
y_model_ic <- lm(pseudo_ic ~ Msim0 + Csim0)

#estimate E[Y^(0, I = 1)] ------------------------------------------------

ey0i_ic <- predict(y_model_ic, 
                   newdata = data.frame(Msim0 = Msim[Asim == 1], 
                                        Csim0 = Csim[Asim == 1]))


#indirect and direct effects
ie_ic <- mean(ey0i_ic) - theta0_ic
de_ic <- theta1_ic - mean(ey0i_ic) 


# Total effect midpoint ----------------------------------------------------------------

sfit_strat_mid <- survfit(Surv(tdata$tmid, 
                               event = tdata$event_mid) ~ tdata$Asim)

rmsts_mid <- summary(sfit_strat_mid, rmean = tstar)$table

theta0_mid <- rmsts_mid[1 ,"rmean"]
theta1_mid <- rmsts_mid[2 ,"rmean"]

total_mid <- theta1_mid - theta0_mid

#outcome model 

sfit_mid <- survfit(Surv(tdata$tmid, 
                         event = tdata$event_mid) ~ 1)

tdata$psvals_mid <- pseudo(sfit_mid, 
                           times = tstar, 
                           type = "rmst")

#fit model on pseduo-observations

y_model_mid <- lm(psvals_mid ~ Msim + Csim,
                  subset = which(Asim == 0),
                  data = tdata)

# estimate E[Y^(0, I = 1)]

ey0i_mid <- predict(y_model_mid, 
                    newdata = tdata[Asim == 1,])

#indirect and direct effects

ie_mid <- mean(ey0i_mid) - theta0_mid
de_mid <- theta1_mid - mean(ey0i_mid)

# All results -------------------------------------------------------------

result <- c(nsample, 
            total_rc, ie_rc, de_rc,
            total_ic, ie_ic, de_ic,
            total_mid, ie_mid, de_mid)

