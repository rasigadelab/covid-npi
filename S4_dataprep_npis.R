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
# NON-PHARMACEUTICAL INTERVENTIONS

source("helpers.R")

load("rdata/country_baseline_200512.Rdata")
load("rdata/branchlengths_200512.Rdata")

##################################################################
# SELECT BRANCHES WITH COUNTRY LIKELIHOOD >95%

blct <- bl[ country_lik >= 0.95 ]

# Add delays relative to root and relative to first divergence in country
blct[, first_div := min(time_source), by = country]
blct[, delay_global := (time_source - min(time_source))*364.25]
blct[, delay_in_country := (time_source - first_div)*364.25]

##################################################################
# INTERVENTIONS DATA

# Sanitize data
interventions <- fread("data/OxCGRT_latest_200512.csv")
names(interventions)

interventions[, Date := ymd(Date)]

# Retain indicators involved in the stringency response index, ie C1-8 and H1

itvnames <- matrix(c(
  "CountryName", "country",
  "Date", "itv_date",
  "C1_School closing", "itv_school",
  "C2_Workplace closing", "itv_workplace",
  "C3_Cancel public events", "itv_public_events",
  "C4_Restrictions on gatherings", "itv_gatherings",
  "C5_Close public transport", "itv_public_transport",
  "C6_Stay at home requirements", "itv_stay_at_home",
  "C7_Restrictions on internal movement", "itv_movement",
  "C8_International travel controls", "itv_intl_travel",
  "H1_Public information campaigns", "itv_info_campaigns",
  "StringencyIndex", "itv_stringency"
), ncol = 2, byrow = T)

setnames(interventions, itvnames[,1], itvnames[,2])
interventions <- interventions[, itvnames[,2], with = F]

# Align country names
setdiff(ctpred$country, interventions$country)

repairs <- matrix(c(
  # "", "Cambodia",
  "Democratic Republic of Congo", "Democratic Republic of the Congo",
  # "", "Georgia",
  # "", "Latvia",
  # "", "Nepal",
  # "", "Senegal",
  "Slovak Republic", "Slovakia",
  "Russia", "Russian Federation",
  "Brunei", "Brunei Darussalam"
), ncol = 2, byrow = T)

unique(interventions[grep("Lat", country), .(country)])

apply(repairs, 1, function(r) {
  interventions[ country == r[1], country := r[2]]
})

# Discard countries not present in dataset
interventions <- interventions[country %in% unique(blct$country)]

names(interventions)

# Examine ordinal-scale interventions if needed, before binarizing
if(FALSE) {
  itv_var <- "itv_info_campaigns"
  dcast_formula <- formula(sprintf("country ~ %s", itv_var))
  itv_by_country <- dcast( interventions[ !is.na(interventions[[itv_var]]) & interventions[[itv_var]] != 0 ], dcast_formula, fun.aggregate = function(x) (length(x) > 0) + 0, value.var = "country")
  colSums(itv_by_country[, -c("country")])
  
  table(interventions[ !is.na(itv_school), .(country, itv_school)])
}


# Binarize interventions
interventions[ , itv_school := (itv_school >= 2) + 0]
interventions[ , itv_workplace := (itv_workplace >= 2) + 0]
interventions[ , itv_public_events := (itv_public_events >= 2) + 0]
interventions[ , itv_gatherings := (itv_gatherings >= 2) + 0]
interventions[ , itv_public_transport := (itv_public_transport >= 2) + 0]
interventions[ , itv_stay_at_home := (itv_stay_at_home >= 2) + 0]
interventions[ , itv_movement := (itv_movement >= 2) + 0]
interventions[ , itv_intl_travel := (itv_intl_travel >= 2) + 0]
interventions[ , itv_info_campaigns := (itv_info_campaigns >= 2) + 0]

# Save ordered variable and label names for consistency
itv_vars <- c("itv_info_campaigns", "itv_intl_travel", "itv_school",
              "itv_public_events", "itv_gatherings", "itv_workplace",
              "itv_movement", "itv_public_transport", "itv_stay_at_home")
itv_labels <- c("Information campaign", "Restrict intl. travel", "Education lockdown",
                "Cancel public events", "Restrict gatherings >100 pers.", "Close workplaces", 
                "Restrict internal movements", "Close public transport", "Home containment")
cbind(itv_vars, itv_labels)
itv_varlist <- as.list(itv_vars)
names(itv_varlist) <- itv_vars

##################################################################
# EXPORT

save(blct, interventions, itv_vars, itv_labels, itv_varlist, file = "rdata/interventions_200518.Rdata")
