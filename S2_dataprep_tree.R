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
# TREE PREPARATION SCRIPT

source("helpers.R")
load("rdata/country_baseline_200512.Rdata")

# Re-run computation-hungry models? (SLOW)
rebuild <- FALSE

########################################################################################################
# LOAD NEXTSTRAIN TREE AND METADATA

# Metadata contain country and sampling date.
# The unit of time in the Nextstrain tree is a year. Use aligndate/dealigndate
# in 'helpers.R' to convert from/to numeric representation
meta_raw <- fread("data/nextstrain_ncov_global_selected_metadata_5211_200512.tsv")
tree_raw <- read.tree("data/nextstrain_ncov_global_timetree_5211_200512.nwk")

plot(tree_raw, show.tip.label = F)


########################################################################################################
# Align metadata with tree labels

meta <- data.table(meta_raw)
tree <- tree_raw

stopifnot(length(setdiff(tree_raw$tip.label, meta_raw$Strain)) == 0)
stopifnot(length(setdiff(meta_raw$Strain, tree_raw$tip.label)) == 0)

# Order meta with tree
meta <- data.table(merge(data.table(Strain = tree$tip.label), meta, by = "Strain", sort = F))
setnames(meta, "Strain", "strain")
stopifnot(all(meta$strain == tree$tip.label))

########################################################################################################
# ALIGN COUNTRY NAMES WITH COUNTRY_BASELINE STANDARD

setnames(meta, "Country", "country")
table(meta$country, useNA = "a")
length(unique(meta$country))

# Align country names with ISO3166
setdiff(meta$country, country_baseline$country)
country_baseline[grepl("Brunei", country)]

# Need to average out intra-country values
repairs <- matrix(c(
  "Russia", "Russian Federation",
  "USA", "United States",
  "Brunei", "Brunei Darussalam"
), ncol = 2, byrow = T)

apply(repairs, 1, function(r) {
  meta[country == r[1], country := r[2]]
})

stopifnot(length(setdiff(meta$country, country_baseline$country)) == 0)

########################################################################################################
# EXCLUDE NON-HUMAN SAMPLES

# Tree can contain non-human samples
table(meta$Host, useNA = "a")

im_human <- meta$Host == "Human"

meta <- meta[im_human]
root.time <- tree$root.edge
tree <- drop.tip(tree, which(!im_human))
tree$root.edge <- root.time

tree
stopifnot(all(meta$strain == tree$tip.label))

###########################################################
# EXPORT

save(meta_raw, tree_raw, meta, tree, file = "rdata/nextstrain_tree_200512.Rdata")