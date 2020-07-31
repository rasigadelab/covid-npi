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

########################################################################################################
#' HELPER ROUTINES
#' 

library(data.table)
library(ape)
library(lubridate)
library(survival)
library(coxme)
library(stringr)
library(beeswarm)
library(plotrix)
library(wordcloud)

# DATE MANAGEMENT

#' Transform date to numeric compatible with GISAID year-based coding
aligndate <- function(d) as.numeric(2020 + (d - ymd("2020-01-01")) / 364.25)

#' Recode aligned, numeric date to standard date
dealigndate <- function(x) {
  ymd("2020-01-01") + (x - 2020)*364.25
}

# I/O MANAGEMENT

#' Automatic file name extension for SVG file
svgf <- function(file, width, height)   eval( {
  svg(sprintf("%s.svg", file), width, height)
}, envir = .GlobalEnv)

#' Automatic file name extension for SVG file
pdff <- function(file, width, height)   eval( {
  cairo_pdf(sprintf("%s.pdf", file), width, height)
}, envir = .GlobalEnv)

#' Export to Excel via the clipboard
xlclipboard <- function(x, dec = ".") {
  if(is.matrix(x)) {
    x <- data.table(row = rownames(x), x)
  }
  write.table(x, file = "clipboard", sep = "\t", na = "", row.names = FALSE, dec = dec)
}

# SEQUENCE MANAGEMENT

#' Fill NA values in an ordered sequence, starting with 0 and repeating the last
#' valid value for other NAs
repair_na_sequence <- function(x, start_value = 0) {
  # Repair first NAs in series
  first_non_na <- min(which(!is.na(x)))
  if(!is.finite(first_non_na)) return(rep(start_value, length(x)))
  if(first_non_na > 1) x[1:(first_non_na - 1)] <- start_value
  # Repair next NAs in series
  while(any(is.na(x))) {
    first_na <- min(which(is.na(x)))
    # Check for end of sequence
    if(all(is.na(x[first_na:length(x)]))) {
      length_na <- length(x) - first_na + 1
    } else {
      length_na <- min(which(!is.na(x[first_na:length(x)]))) - 1
    }
    x[first_na:(first_na + length_na - 1)] <- x[first_na - 1]
  }
  x
}

# PHYLOGENETIC TREE MANAGEMENT

#' Force tree to binary with non-zero branches
#' See http://blog.phytools.org/2018/05/when-phylogeny-fails-isbinary-but-is.html
tree_force_binary <- function(tree, eps = 1e-5) {
  treeb <- multi2di(tree, random = F)
  treec <- collapse.singles(treeb, root.edge = T)
  treec$edge.length[treec$edge.length < eps] <- eps
  stopifnot(is.binary(treec))
  return(treec)
}

# SURVIVAL PLOT

# Make figures with survfit
survplot <- function(svf,
                     stratum = "Intervention",
                     labels = c("No", "Yes"),
                     line.col = c(rgb(0.6, 0, 0), rgb(0,0,0.6)),
                     confint.col = c(rgb(0.6,0,0,0.2), rgb(0,0,0.6,0.2)),
                     legend.pos = "bottomleft",
                     xlim) {
  svfs <- summary(svf)
  stratalev <- levels( svfs$strata )
  
  if(missing(xlim)) {
    xl <- range(svf$time)
  } else {
    xl <- xlim
  }
  
  # xl <- if(missing(xlim) range(svf$time) else xlim
  plot(0, type = "n", xlim = xl, ylim = c(0, 1),
       xlab = "Time (days)\n",
       ylab = "Proportion without\ntransmission", las = 1)
  for(i in 1:length(stratalev)) {
    par(lend = 1)
    # ft <- (svfs$strata == stratalev[i])
    ft <- (1:length(svf$time)) <= cumsum(svf$strata)[i]
    if(i > 1) ft <- ft & (1:length(svf$time)) > cumsum(svf$strata)[i - 1]
    t <- svf$time[ ft ]
    y <- svf$surv[ ft ]
    ylo <- svf$lower[ ft ]
    yup <- svf$upper[ ft ]
    lines(c(0, t), c(1, y), lwd = 1, type = "s", col = line.col[i], lend = 1)
    # Confidence interval (use rep(..., each = 2) trick for stepped lines)
    polx <- rep(c(0, t, rev(t)), each = 2)[-1]#[-(4*sum(ft) + 2)]
    poly <- rep(c(1, ylo, rev(yup[-length(yup)]), 1), each = 2)[-(4*sum(ft) + 2)]
    polygon(polx, poly, lwd = 0.5, border = NA, col = confint.col[i], lend = 1)
    # Censor marks
    ftcensor <- svf$n.censor[ ft ] > 0
    points(t[ ftcensor ], y[ ftcensor ], pch = 3, col = line.col[i])
  }
  # Count events and records
  e0 <- svfs$table[1,4]
  n0 <- svfs$table[1,3]
  e1 <- svfs$table[2,4]
  n1 <- svfs$table[2,3]
  # fstr <- "%s, %g / %g +"
  fstr <- "%s, n=%g (%g+)"
  legend(legend.pos, bty = "n", title = stratum, legend = c(
    sprintf(fstr, labels[1], e0, n0 - e0),
    sprintf(fstr, labels[2], e1, n1 - e1)
  ), fill = confint.col)
}

# SURVIVAL MODEL MANAGEMENT

#' Extract fixed coefficients from coxme.object in a format compatible
#' with summary(coxph.object)
coxme_coef <- function(x) {
  beta <- x$coefficients
  nvar <- length(beta)
  nfrail <- nrow(x$var) - nvar
  se <- sqrt(diag(x$var)[nfrail + 1:nvar])
  tmp <- cbind(beta, exp(beta), se, beta/se,
               # 1 - pchisq((beta/se)^2, 1)
               2*exp(pnorm(-abs(beta/se), log = T))
  )
  dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", 
                                       "se(coef)", "z", "p"))
  tmp
}

#' Express HR in percentage change
to_percents <- function(coefs) {
  t <- data.table(coef = as.double(coefs[[1]]))
  t[, coef := (exp(coef) - 1) * 100]
  t[, lower := (exp(coefs[[1]] - 1.96 * coefs[[3]]) - 1) * 100]
  t[, upper := (exp(coefs[[1]] + 1.96 * coefs[[3]]) - 1) * 100]
  return(t[])
}

#' Express HR with 95%CI
to_hr <- function(coefs) {
  t <- data.table(coef = as.double(coefs[[2]]))
  t[, coef := exp(coef)]
  t[, lower := exp(coefs[[1]] - 1.96 * coefs[[3]])]
  t[, upper := exp(coefs[[1]] + 1.96 * coefs[[3]])]
  return(t[])
}

#' Base plot colors
stdcols <- list(
  dark = c(
    rgb(0.6, 0, 0),
    rgb(0, 0, 0.6),
    rgb(0.8, 0.6, 0),
    rgb(0,0.6,0)
  ),
  light = c(
    rgb(0.9, 0.6, 0.6),
    rgb(0.6,0.6,0.9),
    rgb(1,0.8,0.4),
    rgb(0.6,0.9,0.6)
  ),
  transp = c(
    rgb(0.6, 0, 0, 0.3),
    rgb(0, 0, 0.6, 0.3),
    rgb(0.8, 0.6, 0, 0.3),
    rgb(0,0.6,0,0.3)
  )
)