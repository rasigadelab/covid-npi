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
# GENERATE SUPPLEMENTARY TABLE WITH COUNTRY INFORMATION

source("helpers.R")

print(load("rdata/nextstrain_tree_200512.Rdata"))
print(load("rdata/branchlengths_200512.Rdata"))
print(load("rdata/country_baseline_200512.Rdata"))
print(load("rdata/interventions_200518.Rdata"))

###########################################################################################
# Distribution of sequences across countries
tb <- table(meta$country)
length(tb)

country_table <- data.table(country = names(tb), n_genomes = as.integer(tb))[order(country)]

# Add no. of assigned internal and terminal branches
branchinfo <- blct[, .(n_internal = sum(tip == FALSE), n_terminal = sum(tip == TRUE)), by = country]
country_table <- merge(country_table, branchinfo, by = "country", all.x = T, sort = F)

# Add maximal stringency index
max_stringency <- interventions[, .(max_stringency = max(itv_stringency, na.rm = T)), by = country][is.finite(max_stringency)]
country_table <- merge(country_table, max_stringency, by = "country", all.x = T, sort = F)

###########################################################################################
# Timing of divergence events, declared cases and deaths

# Add date of first divergence
first_edge <- blct[ , .(first_edge = min(first_div)), by = country]
country_table <- merge(country_table, first_edge, by = "country", all.x = T, sort = F)
country_table[, first_edge := dealigndate(first_edge)]

# Merge with dates of 10th case and 10th death
country_table <- merge(country_table, country_baseline[, .(country, date_10cases, date_10deaths)], by = "country", all.x = T, sort = F)

###########################################################################################
# Intervention dates and first edge date (for delay)

itv_onset <- data.table(country = sort(unique(interventions$country)))
for(itv_var in itv_vars) {
  dcast_formula <- formula(sprintf("country ~ %s", itv_var))
  itv_tmp <- dcast(interventions, dcast_formula, value.var = "itv_date",
                   fun.aggregate = function(x) min(x, na.rm = T))[, c("country", "1"), with = F]
  itv_onset <- merge(itv_onset, itv_tmp, by = "country", all.x = T, sort = F)
  setnames(itv_onset, "1", itv_var)
}

country_table <- merge(country_table, itv_onset, by = "country", all.x = T, sort = F)

###########################################################################################
# EXPORT

save(country_table, file = "rdata/country_table_200518.Rdata")

# Export to Excel - Supplementary Table XXX
# Clean row labels
names(country_table)
setnames(country_table, names(country_table), c(
  "Country", "Sample size", "No. internal branches", "No. terminal branches",
  "Max. stringency index", "First divergence", "10th declared case", "10th declared death",
  itv_labels
))
xlclipboard(country_table)
