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
#' 


########################################################################################################
# DATA PREPARATION SCRIPT
# Data alignment and repairs between data sources

source("helpers.R")

########################################################################################################
# ISO3166 LIST OF COUNTRIES

library(maps)
data("world.cities")

# Use 'country' as primary key
ct <- data.table(iso3166)
ct[, country := ISOname]
setkey(ct, 'country')

# DUplicates exist due to map names
ct <- ct[!duplicated(country)]
stopifnot(all(!duplicated(ct$country)))

ct[grepl("Venez", country)]

# Repair important names
repairs <- matrix(c(
  "Republic of Korea", "South Korea",
  "Iran, Islamic Republic of", "Iran",
  "United Kingdom of Great Britain and Northern Ireland", "United Kingdom",
  "Venezuela, Bolivarian Republic of", "Venezuela",
  "Syrian Arab Republic", "Syria"
), ncol = 2, byrow = T)

apply(repairs, 1, function(r) {
  ct[ country == r[1], country := r[2]]
})

##################################################################
# OWID EPIDEMIOLOGY DATA

# stop()

# load OWID data
owid <- fread("data/owid-covid-data_200512.csv", encoding = "UTF-8", na.strings = "")
owid[, date := ymd(date)]
names(owid)

# align country names
setnames(owid, "location", "country")
range(owid$date)

setdiff(owid$country, ct$country)

repairs <- matrix(c(
  "Bolivia", "Bolivia, Plurinational State of",
  "Democratic Republic of Congo", "Democratic Republic of the Congo",
  "Congo", "Republic of Congo",
  "Russia", "Russian Federation"
), ncol = 2, byrow = T)

ct[grep("Venez", country)]

apply(repairs, 1, function(r) {
  owid[ country == r[1], country := r[2]]
})

# Last date per country
owid <- owid[order(-date)]
owid_last <- owid[!duplicated(country)]
range(owid_last$date)

# Prepend names
setnames(owid_last, names(owid_last), paste("owid_", names(owid_last), sep = ""))
names(owid_last)
setnames(owid_last, "owid_country", "country")

# Looks OK, merge
ct <- merge(ct, owid_last, by = "country", all.x = T, sort = F)
stopifnot(all(!duplicated(ct$country)))

###########################################################
# OWID DATES

# Date of 10th death
owid_10death <- owid[ order( country, total_deaths, date ) ][total_deaths >= 10]
owid_10death <- owid_10death[ !duplicated(owid_10death[, .(country)])][order(date), .(country, date_10deaths = date)]

# Date of 10th case
owid_10cases <- owid[ order( country, total_cases, date ) ][total_cases >= 10]
owid_10cases <- owid_10cases[ !duplicated(owid_10cases[, .(country)])][order(date), .(country, date_10cases = date)]

ct <- merge(ct, owid_10death, by = "country", all.x = T, sort = F)
ct <- merge(ct, owid_10cases, by = "country", all.x = T, sort = F)
stopifnot(all(!duplicated(ct$country)))

###########################################################
# OWID CUMULATIVE CASES AND DEATHS
lastcounts <- owid[][order(-date)][!duplicated(country)]

ct <- merge(ct, lastcounts, by = "country", all.x = T, sort = F)

###########################################################
# EXPORT

country_baseline <- ct[, .(country, date_10deaths, date_10cases, total_cases, total_deaths, total_cases_per_million, total_deaths_per_million)]

save(country_baseline, file = "rdata/country_baseline_200512.Rdata")
