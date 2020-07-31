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

#########################################################################
# BUILD TABLE OF BRANCH LENGTHS AND COUNTRIES

source("helpers.R")
load("rdata/country_baseline_200512.Rdata")
load("rdata/nextstrain_tree_200512.Rdata")

# Re-run computation-hungry models? (SLOW)
rebuild <- FALSE

#########################################################################
# ACE ON COUNTRY

#' Perform ACE on forced binary tree with infinitesimal non-zero branch lengths
#' where required. 
#' 

# Branches less than a day ?
aligned_day <- 1/364.25
table(tree$edge.length >= aligned_day)
summary(tree$edge.length[tree$edge.length >= aligned_day])
summary(tree$edge.length[tree$edge.length <  aligned_day])

# Infinitesimal branch length = 1hour
branch_min_length <- aligned_day/24

# Define arbitrary minimal branch length as 1 hour
treed <- tree_force_binary(tree, branch_min_length)
summary(treed$edge.length)

# Ancestral state reconstruction -- SLOW!
if(rebuild) {
  table(meta$country, useNA = "a")
  xstate_full <- ace(meta$country, treed, type = "discrete")
  apply(xstate_full$lik.anc, 2, sum)
  save(xstate_full, file = "rdata/xstate_200512.Rdata")
} else {
  load("rdata/xstate_200512.Rdata")
}

#####################################################################
# BUILD BRANCH LENGTH TABLE

# Build branch length table
n <- Ntip(treed)
bl <- data.table(source = treed$edge[,1], target = treed$edge[,2], length = treed$edge.length)
bl[ , tip := target <= n]
# Reorder by source
setkey(bl, "source")

# Starting point (root time)
root.time <- tree$root.edge

# Include time elapsed with each branch
heights <- node.depth.edgelength(treed)
length(heights)
bl[, time_source := heights[source] + root.time]
bl[, time_target := heights[target] + root.time]
bl[, time_center := 0.5*(time_source + time_target)]

# Initial internal and terminal branches
c(Nedge(tree) - Ntip(tree), Ntip(tree))
# How many branches to discard due to numerically zero lentgh ?
table(bl[ , .(tip, zero = length <= branch_min_length)])

#########################################################################
# ALIGN COUNTRY PROBABILITIES WITH BRANCH LENGTH TABLE

# Align matrix of probabilities, where each row is the edge source,
# hence always an internal node
cty <- xstate_full$lik.anc
cty_mat <- cty[bl$source - n,]
dim(cty_mat)
dim(bl)

colnames(cty_mat)

##################################################################
# ALIGN COUNTRY-WISE PREDICTORS WITH ANCESTRAL STATES

# Table of country_wise predictors aligned with the ancestral likelihoods
ctpred <- merge(data.table(country = colnames(cty_mat)), country_baseline, all.x = T, sort = F)
stopifnot(all(colnames(cty_mat) == ctpred$country))

##################################################################
# ADD BEST-MATCH COUNTRY TO BRANCH LENGTH TABLE

# Distribution of most likely country ?

max_lik <- apply(cty_mat, 1, max)
hist(max_lik)
mean(max_lik > 0.95)

boxplot(bl$time_source ~ max_lik > 0.95)
boxplot(bl$length ~ max_lik > 0.95)

# OK, keep countries as factor and discard unresolved nodes
best_match <- apply(cty_mat, 1, which.max)
best_match <- factor(best_match, levels = 1:ncol(cty_mat))
levels(best_match) <- colnames(cty_mat)

table(best_match)

bl[, country := best_match]
bl[, country_lik := max_lik]

# Define a branch ID that will be used to identify sub-intervals and to control for non-independence in ME models
bl[, id := 1:nrow(bl)]

###########################################################
# EXPORT

save(treed, bl, cty_mat, ctpred, file = "rdata/branchlengths_200512.Rdata")
