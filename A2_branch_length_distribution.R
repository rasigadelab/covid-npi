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

##########################################################################
# BRANCH LENGTH DISTRIBUTION ANALYSIS
 
source("helpers.R")

print(load("rdata/nextstrain_tree_200512.Rdata"))
print(load("rdata/branchlengths_200512.Rdata"))
print(load("rdata/country_baseline_200512.Rdata"))
print(load("rdata/interventions_200518.Rdata"))
print(load("rdata/country_table_200518.Rdata"))

# Re-run computation-hungry code? (SLOW)
rebuild <- FALSE

####################################################################
# TIME SCALES Express length in days and delays in weeks

blct[, length := length * 364.25]
blct[, delay_global := delay_global / 7]
blct[, delay_in_country := delay_in_country / 7]

# BL distribution
table(blct$tip)
summary(blct[tip == FALSE]$length)
summary(blct[tip == TRUE]$length)

# blct[tip == FALSE, .(mean = )]

# Plot tree with colored edges
edgecolor <- c(rgb(0.6,0,0), rgb(0,0,0.6))[ 1 + (treed$edge[,2] <= Ntip(treed)) ]

if(rebuild) {
  pdff("fig/radialtree", 8, 8)
  plot(treed, show.tip.label = FALSE, type = "fan",
       edge.color = edgecolor, open.angle = 10, rotate.tree = 5)
  # add.scale.bar(length = 31/364.25, lwd = 2, labels = "31 days")
  dev.off()
}

if(rebuild) {
  svgf("fig/branchlength_histogram", 4, 4)
  {
    par(mfrow = c(1,1))
    par(mar = c(4,4,4,4))
    breaks <- seq(0, 100, by = 5)
    hist(blct[tip == FALSE]$length, breaks = breaks, border = rgb(0.6,0,0), col = rgb(0.6,0,0,0.3), main = "", xlab = "Branch length (days)")
    hist(blct[tip == TRUE]$length, breaks = breaks, border = rgb(0,0,0.6), col = rgb(0,0,0.6,0.3), add = T)
    legend("topright", bty = "n", title = "Branch type", legend = c("internal (time-to-event)", "terminal (censored)"), 
           fill = c(rgb(0.6,0,0,0.3), rgb(0,0,0.6,0.3)))
  }
  dev.off()
}
