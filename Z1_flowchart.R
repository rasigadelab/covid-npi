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

#####################################################################
# FLOWCHART
# Summarize inclusion/exclusion of data
 
source("helpers.R")

print(load("rdata/nextstrain_tree_200512.Rdata"))
print(load("rdata/branchlengths_200512.Rdata"))
print(load("rdata/country_baseline_200512.Rdata"))
print(load("rdata/interventions_200518.Rdata"))
print(load("rdata/country_table_200518.Rdata"))

#####################################################################
# BRANCH COUNTS

branchcount <- function(tree) {
  c(terminal = Ntip(tree), internal = Nnode(tree) - 1, total = Ntip(tree) + Nnode(tree) - 1)
}

branchcount(tree_raw)
# Exclude non-human samples
branchcount(tree_raw) - branchcount(tree)
branchcount(tree)
# Binarize polytomies
branchcount(treed) - branchcount(tree)
branchcount(treed)

table(bl$tip)
# Exclude poorly resolved countries
diff <- table(bl$tip) - table(blct$tip)
c(diff, total = sum(diff))
tb <- table(blct$tip)
c(tb, total = sum(tb))

# Branches in countries with known NPIs
tb <- table(blct[ country %in% country_table[!is.na(itv_info_campaigns)]$country ]$tip)
diff <- table(blct$tip) - tb
c(diff, total = sum(diff))
c(tb, total = sum(tb))

# Very few branches in countries with missing data, this is because most countries
# were excluded at the previous step

length(unique(blct$country))

#####################################################################
# REPRESENTATIVENESS OF EXCLUDED BRANCHES

bl[, uncertain_country := country_lik < 0.95]
bl[, repr := ""]
bl[
  tip == TRUE & uncertain_country == FALSE, repr := "Terminal\n(included)"
]
bl[
  tip == TRUE & uncertain_country ==  TRUE, repr := "Terminal\n(excluded)"
]
bl[
  tip == FALSE & uncertain_country == FALSE, repr := "Internal\n(included)"
]
bl[
  tip == FALSE & uncertain_country ==  TRUE, repr := "Internal\n(excluded)"
  ]

pdff("fig/excluded_branches", 6, 4)
{
  par(mar = c(4,10,2,2))
  boxplot(length * 364.25 ~ repr, bl, horizontal = TRUE, log = "", xlab = "Branch length (days)\n", ylab = "",
          las = 1, col = rgb(0.6,0,0,0.3))
  mtext("Branch type", side = 3, adj = 0, at = -40)
}
dev.off()

#####################################################################
# COUNTRY COUNTS

# Countries contributing >1 genome
names(country_table)
nrow(country_table)

# Countries with no assigned branch
table(is.na(country_table$n_internal))

# Countries with >= 10 internal branches (Figure 2)
table(country_table$n_internal >= 10)

# Countries with known NPIs
nrow(country_table) - sum(!is.na(country_table$itv_info_campaigns))
sum(!is.na(country_table$itv_info_campaigns))





