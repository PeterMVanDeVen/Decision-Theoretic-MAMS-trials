# Supplementary R code for manuscript "Bayesian adaptive decision-theoretic 
# designs for multi-arm multi stage clinical trials" published in 
# Statistical Methods in Medical Research
#
# Authors: A. Bassi, J. Berkhof, D. de Jong & P.M. van de Ven
#
# Correspondence/enquiries: p.vandeven@amsterdamumc.nl



# This script requires that the R script 'Functions for selecting the best from
# three, four or five experimental arms.R' was previously run and packages snow, 
# doSNOW and doParallel for parallel computing were installed.
# The required R script can be downloaded from 
# https://github.com/PeterMVanDeVen/Decision-Theoretic-MAMS-trials



# This example illustrates the use of the following functions for the 
# setting with five experimental arms: 
#
# create.d  : Simulates trial data
# prdrop    : Evaluates trials with adaptive stopping when arms can be dropped
# prnodrop  : Evaluates trials with adaptive stopping without dropping of arms
# reevaluate: Reevaluates trials using output from prdrop and prnodrop
#             for higher thresholds gamma (=C/Q) 
# pcalc     : Calculates operating characteristics (samples sizes and final decisions) 
#             using output from prdrop, prnodrop and reevaluate
#
# The required functions provided in a separate R script 'Functions for selecting the best from
# three, four or five experimental arms.R'



# Packages snow, doSNOW an doParallel are required for parallel 
# computing 

library(snow)
library(doSNOW)
library(doParallel)

# Specify the number of trials to be simulated (simulation size)

Ntrials <- 1000

# Specify the response rate vector under which the frequentist 
# properties of the trial are to be evaluated

resp <- c(0.2,0.7,0.8,0.2,0.7)

# Specify the cap (= maximum number of subjects included in the trial)
# In this example the cap is set at a high number that is not reached 
# for the response rate vector and threshold (C/Q) used for stopping

maxpt <- 1000

# Set margin delta at 0 for symmetric setting with only experimental
# arms

delta <- 0

# Set dropctrl at TRUE to allow dropping of all experimental arms 

dropctrl <- TRUE

# Specify batch size for stage 1 
# 'burn' should be specified as the total number of subjects in the stage
# divided by the number of arms at start of the trial
# 'burn <- 4' here corresponds to 20 subjects in the first stage

burn <- 4

# Specify the batch size for stage 2, 3, etc. 
# 'batch' should be specified as the total number of subjects per stage
# divided by the number of arms at start of the trial
# 'batch <- 4' here corresponds to 20 subjects in each of the stages 
# 2, 3, ...

batch <- 4

# Specify the minimally required increase in probability of a correct
# decision that is required for continuation
# Note that in the manuscript this threshold gamma is referred to as 
# C/Q

gamma <- 1/1000

# Set the parameters a and b for the Beta(a,b) prior 
# The prior is common for all arms
# We use a uniform prior and set a = b = 1

prior <- c(1,1)

# Generate data for (Ntrials) trials using the create.d function

Ydata <- create.d(resp, maxpt, Ntrials)

# Evaluate all trials simulated using the prdrop function for the
# setting where dropping of arms is allowed
# This is the most computationally intensive part. Under settings considered and using 8 cores on a i5
# processor laptop a call to prdrop required approximely 4 hours computation time (1000 trials)
# Computation increases with smaller stages sizes and smaller gamma

trial.out <- prdrop(Ydata, delta, gamma, burn, Ntrials, batch, prior, maxpt, dropctrl)

# Evaluate all trials simulated using the prnodrop function for the
# setting where dropping of arms is not allowed
# Under settings considered and using 8 cores on a i5
# processor laptop a call to prnodrop required approximely 1 hour and 30 minutes 
# computation time (1000 trials)

trialnodrop.out <- prnodrop(Ydata, delta, gamma, burn, Ntrials, batch, prior, maxpt)

# Summarize the results using the pcalc function 
# Function returns the probability of each final decision and a 
# vector containing total sample sizes for all trials

res_drop   <- pcalc(trial.out, prior, delta, nr_dec = 5)
res_nodrop <- pcalc(trialnodrop.out, prior, delta, nr_dec = 5)

# Look at the results

res_drop
res_nodrop

# Look at proportion of each final decision
# res_drop$dec and res_nodrop$dec are vectors (p1,p2,p3,p4,p5) where
# p1: proportion of trials in which experimental arm 1 is selected
# p2: proportion of trials in which experimental arm 2 is selected
# p3: proportion of trials in which experimental arm 3 is selected 
# p4: proportion of trials in which experimental arm 4 is selected 
# p5: proportion of trials in which experimental arm 5 is selected 


res_drop$dec
res_nodrop$dec

# Specifically look at proportion of a trials with a correct decision
# (p3)

res_drop$dec[3]
res_nodrop$dec[3]

# res_drop$N and res_nodrop$N contain the trial sizes
# Look at average trial size

mean(res_drop$N)
mean(res_nodrop$N)

# We now show how to reevaluate trials using higher threshold 
# C/Q=100 yielding shorter trials

new_gamma <- 1/500

# We reevaluate the simulate trials with and without dropping with 
# the new threshold C/Q = 1/500

trial.out.new.gamma <- reevaluate(trial.out, Ydata,  new_gamma)
trialnodrop.out.new.gamma <- reevaluate(trialnodrop.out, Ydata, new_gamma)

# We summarize the results for the higher threshold using the 
# pcalc function 

res_drop.new.gamma   <- pcalc(trial.out.new.gamma, prior, delta, nr_dec = 5)
res_nodrop.new.gamma <- pcalc(trialnodrop.out.new.gamma, prior, delta, nr_dec = 5)

# Look at the results for the higher threshold

res_drop.new.gamma 
res_nodrop.new.gamma

# Specifically look at proportion of a trials with a correct decision 
# (p3)

res_drop.new.gamma$dec[3]
res_nodrop.new.gamma$dec[3]



