########################################################################################################
#' Impact of non-pharmaceutical interventions on the spread of SARS-CoV-2:
#' a cross-national, multivariate genomic analysis
#' 
#' Jean-Philippe Rasigade, M.D., Ph.D., Anaïs Barray, M.Sc., Yoann Vigouroux, M.Sc.,
#' Charlène Coquisart, M.Sc., Julie T. Shapiro, Ph.D., Antonin Bal, Pharm.D., Grégory Destrat, Pharm.D.,
#' Philippe Vanhems, M.D., Ph.D., Bruno Lina, M.D., Ph.D., Laurence Josset, Pharm.D, Ph.D,Thierry Wirth, Ph.D.
#' 
#' (c) 2020 Jean-Philippe Rasigade, <jean-philippe.rasigade(at)univ-lyon1.fr>
#' 
#' MIT LICENSE

#########################################################
# EPIDEMIOLOGICAL SIR MODELS

library(deSolve)

source("helpers.R")
print(load("rdata/interventions_200518.Rdata"))

#####################################################
# Load results from full multivariate model
load("rdata/cxme_full_id_200716.Rdata")
coefmat <- coxme_coef( cxme_full_id )

# Compute median intervention delay relative to divergence
load("rdata/delay_table_200518.Rdata")
itv_timing <- delay_table[, .(mean = mean(delay), median = median(delay), q1 = quantile(delay, 0.25), q3 = quantile(delay, 0.75)), by = itv]

# Update intervention labels
itv_labels[5] <- "Restrict gatherings"
itv_timing$itv <- itv_labels

######################################################
# SIR MODEL ROUTINES

# Logistic growth model to illustrate cumulative HR results
loggrowth <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    beta <- ifelse(Time < time_itv, beta, beta_itv)
    
    dS <- - beta * I * S / N
    dI <- beta * I * S / N - gamma * I
    dR <- gamma * I
    # Cumulative cases
    dCum <- beta * I * S / N
    
    return(list(
      c(dS, dI, dR, dCum)
    ))
    
  })
}

# SIR result with intervention
getSIR <- function(timing = 0, coef = 0, length = 4, baseline) {
  suppressWarnings(with(baseline, {
    pars <- c(
      beta = beta,
      gamma = gamma,
      N = N0,
      time_itv = timing / 30, # BEWARE that intervention timing is in days while SIR timing is in months
      beta_itv = beta * exp(coef) # Modified infection rate after intervention
    )
    
    yini <- c(
      S = N0 - I0,
      I = I0,
      R = 0,
      Cum = I0
    )
    
    times <- seq(0, length, by = 0.05)
    # print(pars)
    out_itv <- data.table(ode(yini, times, loggrowth, pars))
    
    # Filter pre-intervention data
    out_itv <- out_itv[time >= timing / 30]
  }))
}

#####################################################
# Baseline parameters (in *months*)
R0 = 3
gamma = 2
beta = R0 * gamma
N0 = 1e6
I0 = 100 # First divergence event

baseline <- list(
  R0 = R0,
  gamma = gamma,
  beta = beta,
  N0 = N0,
  I0 = I0 # First divergence event
)

#####################################################
# Combine NPI timing and coefficients

# Combine with coeffs and CIs
itv_timing[, coef := coefmat[,1]]
itv_timing[, coef_lo := coef - 1.96*coefmat[,3]]
itv_timing[, coef_hi := coef + 1.96*coefmat[,3]]
itv_timing[, coef_se := coefmat[,3]]

sir_len <- 8 # Simulation length in months

sir_ref <- getSIR(0, 0, sir_len, baseline)

svgf("fig/SIR_interventions", 6, 6)
{
  options(scipen = 999) # Penalize using scientific notation in axes
  par(mfrow = c(3,3))
  par(mar = c(5, 5, 2, 2))
  par(xpd = FALSE)
  
  scaling <- 1
  
  for(i in 1:nrow(itv_timing)) {
    sir_med <- getSIR(itv_timing$median[i], itv_timing$coef[i], sir_len, baseline)
    sir_lo <- getSIR(itv_timing$median[i], itv_timing$coef_lo[i], sir_len, baseline)
    sir_hi <- getSIR(itv_timing$median[i], itv_timing$coef_hi[i], sir_len, baseline)
    
    {
      plot(sir_ref$time, sir_ref$Cum * scaling, type = "l", log = "y", ylim = c(I0, N0), las = 1,
           xlab = "Time (months)", lwd = 2, col = "darkgrey", ylab = "Cumulative cases\n")
      lines(sir_med$time, sir_med$Cum * scaling, lwd = 2, col = stdcols$dark[1])
      polygon(
        c(sir_lo$time, rev(sir_lo$time)), c(sir_lo$Cum * scaling, rev(sir_hi$Cum * scaling)),
        col = stdcols$transp[1], lwd = 0.5)
      # Allow title to overflow figure region in multiplot
      par(xpd = NA)
      title(main = itv_timing$itv[i])
      par(xpd = FALSE)
    }
  }
}
dev.off()


##########################################################################
# PLOTS OF SIMULTANEOUS CASES

quantile_lines <- seq(0.025, 0.975, by = 0.025)
quantile_col <- rgb(0.6, 0, 0, 0.075)

pdff("fig/fig_s9_SIR_NPI_daily_cases", 6, 4.5)
{
  options(scipen = 999) # Penalize using scientific notation in axes
  par(mfrow = c(3,3))
  par(mar = c(3.2, 5, 2, 1))
  par(xpd = FALSE)
  par(font.main = 1)
  
  scaling <- 1e-3
  sir_len <- 8 # Simulation length in months
  sir_ref <- getSIR(0, 0, sir_len, baseline)
  
  for(i in 1:nrow(itv_timing)) {
    print(itv_timing$itv[i])
    plot(sir_ref$time, sir_ref$I * scaling, type = "l", log = "", ylim = c(0, 3.2e5 * scaling), las = 1,
         xlab = "Time (months)\n", lwd = 2, col = "black", ylab = "Simultaneous\ncases x 1,000", lty = 1)
    sir_med <- getSIR(itv_timing$median[i], itv_timing$coef[i], sir_len, baseline)
    lines(sir_med$time, sir_med$I * scaling, lwd = 2, col = stdcols$dark[1])
    
    # Quantile lines
    for(qt in quantile_lines) {
      qcoef <- qnorm(qt, itv_timing$coef[i], itv_timing$coef_se[i])
      sir <- getSIR(itv_timing$median[i], qcoef, sir_len, baseline)
      lines(sir$time, sir$I * scaling, lwd = 3, col = quantile_col)
    }
    
    # Allow title to overflow figure region in multiplot
    par(xpd = NA)
    title(main = itv_timing$itv[i])
    par(xpd = FALSE)
  }
}
dev.off()

# Estimate peak amplitude and delay
sir_summary <- function(sir) {
  nmax <- which.max(sir$I)
  sir[nmax]$time
  return(c(
    peak_cases = sir[nmax]$I,
    peak_time = sir[nmax]$time
  ))
}

# Add to main table
itv_timing

cf <- "coef"
timing_col <- "median"

# Simulate longer time for peak estimation
sir_len <- 24

for(cf in c("coef", "coef_lo", "coef_hi")) {
  for(timing_col in c("median", "q1", "q3")) {
    coln <- paste(cf, timing_col, sep = "_")
    print(coln)
    
    res <- data.table(t(
      apply(as.matrix(itv_timing[, c(timing_col, cf), with = F]), 1, function(x) {
        x <- as.double(x)
        sir_summary(getSIR(x[1], x[2], sir_len, baseline))
      })
    ))
    
    setnames(res, names(res), paste(names(res), coln, sep = "_"))
    
    itv_timing <- cbind(itv_timing, res)
  }
}

itv_timing

# Number formatting for large integers
largeint <- function(x) formatC(x, big.mark = ",", format = "f", digits = 0)

# Build proper table

itv_sir_table <- itv_timing[, .(itv)]

itv_sir_table[[ "R0" ]] <- sprintf("%.2f (%.2f to %.2f)",
                                   pmax(0, R0 * exp(itv_timing$coef)),
                                   pmax(0, R0 * exp(itv_timing$coef_lo)),
                                   pmax(0, R0 * exp(itv_timing$coef_hi))
)

for(cf in c("coef", "coef_lo", "coef_hi")) {
  for(timing_col in c("q1", "median", "q3")) {
    
    scaling <- 1000
    cases <- floor(itv_timing[[ paste("peak_cases_coef", timing_col, sep = "_") ]] / scaling)
    cases_lo <- floor(itv_timing[[ paste("peak_cases_coef_lo", timing_col, sep = "_") ]] / scaling)
    cases_hi <- floor(itv_timing[[ paste("peak_cases_coef_hi", timing_col, sep = "_") ]] / scaling)
    
    itv_sir_table[[ paste("peak_cases", timing_col, sep = "_") ]] <-
      sprintf("%s (%s to %s)", largeint(cases), largeint(cases_lo), largeint(cases_hi) )
    
  }
}

for(cf in c("coef", "coef_lo", "coef_hi")) {
  for(timing_col in c("q1", "median", "q3")) {
    
    delay <- itv_timing[[ paste("peak_time_coef", timing_col, sep = "_") ]]
    delay_lo <- itv_timing[[ paste("peak_time_coef_lo", timing_col, sep = "_") ]]
    delay_hi <- itv_timing[[ paste("peak_time_coef_hi", timing_col, sep = "_") ]]
    
    itv_sir_table[[ paste("peak_time", timing_col, sep = "_") ]] <-
      sprintf("%.1f (%.1f to %.1f)", delay, pmin(delay_lo, delay_hi), pmax(delay_lo, delay_hi) )
    
  }
}

itv_sir_table

sir_summary(sir_ref)

xlclipboard(itv_sir_table)


########################################################
# Table S2 - Combinations
# Extract coefficient correlation matrix to estimate standard errors
# Variance of sum is the sum of variances and covariances

m <- fixef(cxme_full_id)
v <- vcov(cxme_full_id)

names(m)
sumrange <- c(1:9)

# Combination function
f <- function(sumrange, perc = TRUE) {
  csum <- sum(m[sumrange])
  se <- sqrt(sum(v[sumrange, sumrange]))
  # Confint of combination
  cfint <- c(coef = csum, lo = csum - 1.96 * se, hi = csum + 1.96 * se)
  
  if(perc == TRUE) cfint <- 100*(exp(cfint) - 1)
  
  return(c(
    # cfint,
    cfint,
    # Stop probability for R0 = 1.5
    stop15 = pnorm(0, log(1.5) + csum, se),
    # Stop probability for R0 = 2.0
    stop20 = pnorm(0, log(2.0) + csum, se),
    # Stop probability for R0 = 2.5
    stop25 = pnorm(0, log(2.5) + csum, se),
    # Stop probability for R0 = 3.0
    stop30 = pnorm(0, log(3.0) + csum, se),
    # Stop probability for R0 = 3.5
    stop35 = pnorm(0, log(3.5) + csum, se)
  ))
}

options(scipen = 1)
t(sapply(1:9, function(n) f(1:n)))

########################################################################
# Implicit restrictions of home containment
# Compute the marginal effect of home containment + implicit restrictions

itv_labels
f(c(4,5,7,9), perc = TRUE)
lockdown_coef <- f(c(4,5,7,9), perc = FALSE)

sapply(1:3, function(i) {
  sir_lockdown <- getSIR(itv_timing$median[9], as.double(lockdown_coef[i]), 12, baseline)
  max(sir_lockdown$I)
})

########################################################################
# Table 2 - Probabilities of reducing Rt below 1 using combined interventions
{
  stoptable <- data.table(t(sapply(1:9, function(n) f(1:n))))
  
  stoptable[, reduction := sprintf("%.1f (%.1f to %.1f)", coef, lo, hi) ]
  stoptable[, stop15_display := sprintf("%.2f", stop15)]
  stoptable[stop15 < 0.01, stop15_display := "<0.01"]
  stoptable[, stop20_display := sprintf("%.2f", stop20)]
  stoptable[stop20 < 0.01, stop20_display := "<0.01"]
  stoptable[, stop25_display := sprintf("%.2f", stop25)]
  stoptable[stop25 < 0.01, stop25_display := "<0.01"]
  stoptable[, stop30_display := sprintf("%.2f", stop30)]
  stoptable[stop30 < 0.01, stop30_display := "<0.01"]
  stoptable[, stop35_display := sprintf("%.2f", stop35)]
  stoptable[stop35 < 0.01, stop35_display := "<0.01"]
}

xlclipboard(stoptable)

########################################################################
# Table S3 - Probabilities of reducing Rt below 1 using individual interventions
{
  stoptable <- data.table(t(sapply(1:9, function(n) f(n))))
  
  stoptable[, reduction := sprintf("%.1f (%.1f to %.1f)", coef, lo, hi) ]
  stoptable[, stop15_display := sprintf("%.2f", stop15)]
  stoptable[stop15 < 0.01, stop15_display := "<0.01"]
  stoptable[, stop20_display := sprintf("%.2f", stop20)]
  stoptable[stop20 < 0.01, stop20_display := "<0.01"]
  stoptable[, stop25_display := sprintf("%.2f", stop25)]
  stoptable[stop25 < 0.01, stop25_display := "<0.01"]
  stoptable[, stop30_display := sprintf("%.2f", stop30)]
  stoptable[stop30 < 0.01, stop30_display := "<0.01"]
  stoptable[, stop35_display := sprintf("%.2f", stop35)]
  stoptable[stop35 < 0.01, stop35_display := "<0.01"]
}

xlclipboard(stoptable)

############################################################
# CUMULATIVE SIR MODELS
# We need different ODE functions to use several change points

# Logistic growth model to illustrate cumulative HR results
loggrowth_cp <- function(Time, State, Pars, breakpoints, coefs) {
  with(as.list(c(State, Pars)), {

    beta <- beta
    # Find smaller change point
    if(Time > breakpoints[1]) {
      cp <- min(which(c(breakpoints, Inf) > Time)) - 1
      
      # print(cp)
      beta <- coefs[cp]
    }

    
    dS <- - beta * I * S / N
    dI <- beta * I * S / N - gamma * I
    dR <- gamma * I
    # Cumulative cases
    dCum <- beta * I * S / N
    
    # if(Time > breakpoints[5]) stop()
    
    return(list(
      c(dS, dI, dR, dCum)
    ))
    
  })
}

# SIR result with several interventions and change points
getSIR_cp <- function(changepoints, coefs, length = 4, baseline) {
  suppressWarnings(with(baseline, {
    pars <- c(
      beta = beta,
      gamma = gamma,
      N = N0
    )
    
    breakpoints <- changepoints / 30 # BEWARE that intervention timing is in days while SIR timing is in months
    coefs <- beta * exp(coefs) # Modified infection rate after intervention
    
    yini <- c(
      S = N0 - I0,
      I = I0,
      R = 0,
      Cum = I0
    )
    
    times <- seq(0, length, by = 0.05)
    
    # ODE function
    fode <- function(Time, State, Pars) loggrowth_cp(Time, State, Pars, breakpoints, coefs)
    
    # print(pars)
    out_itv <- data.table(ode(yini, times, fode, pars))
    
    # Filter pre-intervention data
    out_itv <- out_itv[time >= changepoints[1] / 30]
  }))
}

sir <- getSIR_cp(itv_timing$median, cumsum(itv_timing$coef), sir_len, baseline)
max(sir$I)

# Combination function for coefficients
f <- function(sumrange) {
  csum <- sum(m[sumrange])
  se <- sqrt(sum(v[sumrange, sumrange]))
  # Confint of combination
  cfint <- c(coef = csum, lo = csum - 1.96 * se, hi = csum + 1.96 * se)
  
  return(cfint)
}

cumcoefs <- t(sapply(1:9, function(n) f(1:n)))


############################################################
# CUMULATIVE SIR MODELS - Figure 3

ramp <- colorRampPalette(c("darkgreen", "orange", "darkred"), alpha = 0.5)(nrow(itv_timing))

svgf("fig/cumulative_npi_SIR_curve", 5, 3.5)
{
  options(scipen = 999) # Penalize using scientific notation in axes
  # par(mfrow = c(3,3))
  par(mar = c(5, 5, 2, 2))
  par(xpd = FALSE)
  par(font.main = 1)
  
  scaling <- 1e-3
  sir_len <- 18
   
  sir_ref <- getSIR_cp(0, 0, sir_len, baseline)
  
  plot(sir_ref$time, sir_ref$I * scaling, type = "l", log = "", ylim = c(0, 3.2e5 * scaling), las = 1,
       xlab = "Time (months)\n", lwd = 2, col = "black", ylab = "Simultaneous\ncases x 1,000", lty = 1)
  
  for(i in 1:nrow(itv_timing)) {
    sir_med <- getSIR_cp(itv_timing$median[1:i], cumcoefs[1:i,1], sir_len, baseline)
    lines(sir_med$time, sir_med$I * scaling, lwd = 2, col = ramp[i])
  }
  
  for(i in 1:nrow(itv_timing)) {
    text(sir_len, 3.2e5 * scaling - i*30, itv_labels[i], adj = 1, col = ramp[i])
  }
}
dev.off()
