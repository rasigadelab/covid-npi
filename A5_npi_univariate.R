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
# SEPARATE COX MODELS OF NON-PHARMACEUTICAL INTERVENTIONS

source("helpers.R")

print(load("rdata/nextstrain_tree_200512.Rdata"))
print(load("rdata/branchlengths_200512.Rdata"))
print(load("rdata/country_baseline_200512.Rdata"))
print(load("rdata/interventions_200518.Rdata"))
print(load("rdata/country_table_200518.Rdata"))

# Re-run computation-hungry code? (SLOW)
rebuild <- FALSE

###################################################################
# EXCLUDE INTERVALS IN COUNTRIES WITH UNKNOWN INTERVENTIONS

excluded_countries <- country_table[is.na(itv_info_campaigns)]$country
blct <- blct[ !(country %in% excluded_countries) ]

##################################################################
# BUILD TVC TABLE FOR SINGLE INTERVENTION
# Loop over interventions to derive models with random intercept per country

itv_var <- "itv_school"

make_tvc <- function(itv_var) {
  # Keep only change points for intervention
  # itv <- interventions[ !is.na(interventions[[itv_var]]) ][order(country, itv_date)]
  itv <- data.table(interventions[order(country, itv_date)])
  # Constant intervention name
  setnames(itv, itv_var, "interv")
  
  itv <- itv[, .(country, itv_date, interv)]
  
  # Repair NAs if any
  itv[, interv := repair_na_sequence(interv), by = country]
  
  itv <- itv[!duplicated(itv[, -c("itv_date")])]
  stopifnot(all(!duplicated(itv[, .(country, itv_date)])))
  
  # Transform date to numeric (align with dating in branch length table)
  itv[, itv_date := aligndate(itv_date)]
  
  # Match branch ID with intervention data
  itv <- merge(itv, blct, by = "country", allow.cartesian = TRUE)
  setkey(itv, "itv_date")
  
  tvc <- tmerge(blct, blct, id = id, status = event(time_target - time_source, !tip))
  tvc <- data.table(
    tmerge(tvc, itv, id = id
           , options = list(tdcstart = 0, na.rm = T),
           interv = tdc(itv_date - time_source, interv)
    )
  )
  
  names(tvc)
  tvc <- tvc[, length := (tstop - tstart)*364.25]
  return(tvc)
}

x <- make_tvc("itv_school")
table(x[, .(interv, status)])

###############################################################################
# Global stringency index regression

if(rebuild) {
  cxm <- coxme(Surv(length, status) ~ interv + (1|id) + (1|country), make_tvc("itv_stringency"))
  save(cxm, file = "rdata/coxme_stringency_200715.Rdata")
}
load("rdata/coxme_stringency_200715.Rdata")
cxm
coxme_coef(cxm)
to_percents(coxme_coef(cxm)*10)

##############################################################################
# Individual regressions

# Model fit takes ~25 minutes, parallelize. One file generated per model

# rebuild <- TRUE

if(rebuild) {
  library(parallel)
  
  cl <- makeCluster(length(itv_vars))
  
  clusterEvalQ(cl, {
    library(data.table)
    library(lubridate)
    library(survival)
    library(coxme)
    library(stringr)
    
    source("helpers.R")
    
    print(load("rdata/nextstrain_tree_200512.Rdata"))
    print(load("rdata/branchlengths_200512.Rdata"))
    print(load("rdata/country_baseline_200512.Rdata"))
    print(load("rdata/interventions_200518.Rdata"))
    print(load("rdata/country_table_200518.Rdata"))
  })
  
  clusterExport(cl, c("make_tvc", "interventions", "blct"))
  
  clusterApply(cl, itv_varlist, function(itv_var) {
    save(itv_var, file = sprintf("rdata/flag_varmodel_%s_200715.Rdata", itv_var))
    
    var_model <- coxme(Surv(length, status) ~ interv + (1|id) + (1|country), make_tvc(itv_var))
    
    # SANITY CHECK
    # var_model <- coxme(Surv(length, status) ~ interv + (1|country), make_tvc(itv_var))
    
    save(var_model, file = sprintf("rdata/cxme_varmodel_%s_200715.Rdata", itv_var))
  })
  
  stopCluster(cl)
}

# Gather all models into a model list

model_list <- list()

for(itv_var in itv_vars) {
  load(sprintf("rdata/cxme_varmodel_%s_200715.Rdata", itv_var))
  model_list[[itv_var]] <- var_model
}

model_coefs <- lapply(model_list, coxme_coef)

n_coef <- 1
model_table <- rbindlist(lapply(model_coefs, function(x) {
  data.table(var = rownames(x)[n_coef], t(x[n_coef,]))
}))

model_table$var <- itv_vars
model_table

# Percentage changes
cbind(itv_vars, to_percents(model_table[,-1]))

# Export for combined plot in figure 3
umod_table <- data.table(model_table)
umod_perc <- data.table(cbind(itv_vars, to_percents(model_table[,-1])))

save(umod_table, umod_perc, file = "rdata/npi_univar.Rdata")

#####################################################################
# KAPLAN MEIER PLOT


itv_id <- 1
Sys.setlocale("LC_TIME", "C") # Dates in C standard display
pdff("fig/kaplan_npi", 8, 6)
{
  par(mfrow = c(3,3))
  par(mar = c(4, 5, 1.5, 1))
  for(itv_id in 1:9) {
    itv_var <- itv_vars[itv_id]
    svf <- survfit(Surv(length, status) ~ interv, make_tvc(itv_var))
    survplot(svf, NULL, c("No", "Yes"), legend.pos = "topright", xlim = c(0, 80))
    mtext(itv_labels[itv_id], side = 3, font = 1)
    # Display country-adjusted HR
    hrtext <- {
      cf <- model_table[var == itv_var]$coef
      se <- model_table[var == itv_var]$`se(coef)`
      hrs <-100*(exp(c(cf, cf - 1.96*se, cf + 1.96*se))-1)
      sprintf("Rt (%% change), %.2f (%.2f - %.2f)", hrs[1], hrs[2], hrs[3])
    }

    text(80, 0.1, sprintf("%s", hrtext), adj = c(1,1), cex = 0.8, font = 1)
  }
}
dev.off()


