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
# MULTIVARIATE COX MODELS OF NON-PHARMACEUTICAL INTERVENTIONS

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
# EXTRACT INTERVENTIONS AND SANITIZE NAs

# Keep only change points for *any* intervention
itv <- interventions[order(country, itv_date), -c("itv_stringency")]
itv <- itv[!duplicated(itv[, -c("itv_date")])]
stopifnot(all(!duplicated(itv[, .(country, itv_date)])))

# reorder columns
itv <- itv[, c("country", "itv_date", itv_vars), with = F]

#' There are NAs in the middle of the dataset. Set NA = 0 if first date in country, of NA = last non-NA value otherwise.
itv[, (itv_vars) := lapply(.SD, repair_na_sequence), by = country, .SDcols = itv_vars]

##############################################################
# TABLE S2 - Detailed intervention timeline
xlclipboard(itv)

##################################################################
# BUILD TVC TABLE

# Match branch ID with intervention data
blitv <- data.table(itv)
blitv <- merge(blitv, blct, by = "country", allow.cartesian = TRUE)

# Transform date to numeric (align with dating in branch length table)
blitv[, itv_date := aligndate(itv_date)]
setkey(blitv, "itv_date")

# WARNING a data.table cannot be passed as argument to tmerge, this must be a data.frame

tvc <- tmerge(blct, blct, id = id, status = event(time_target - time_source, !tip))

tvc <- data.table(
  tmerge(tvc, blitv, id = id
         , options = list(tdcstart = 0, na.rm = T),
         itv_info_campaigns = tdc(itv_date - time_source, itv_info_campaigns),
         itv_intl_travel = tdc(itv_date - time_source, itv_intl_travel),
         itv_school = tdc(itv_date - time_source, itv_school),
         itv_public_events = tdc(itv_date - time_source, itv_public_events),
         itv_gatherings = tdc(itv_date - time_source, itv_gatherings),
         itv_workplace = tdc(itv_date - time_source, itv_workplace),
         itv_movement = tdc(itv_date - time_source, itv_movement),
         itv_public_transport = tdc(itv_date - time_source, itv_public_transport),
         itv_stay_at_home = tdc(itv_date - time_source, itv_stay_at_home)
  )
)

names(tvc)
tvc <- tvc[, length := (tstop - tstart)*364.25]


# Distribution of combinations of interventions in subintervals
# Add cumulative and mean length, no. of events
xlclipboard(tvc[ , .N, by = itv_vars][order(-N)])

# Binary distance (proportion of discordant bits) across subintervals
cr <- as.matrix(dist(t(as.matrix(tvc[, (itv_vars), with = F])), method = "binary"))
diag(cr) <- NA

xlclipboard(cr)

# Kendall correlation
kmat <- as.matrix(tvc[, (itv_vars), with = F]) == 1
kmat

kcross <- crossprod(kmat) # Both active
kcrossnot <- crossprod(!kmat) # Neither active
nrow(kmat) - kcross - kcrossnot # XOR

crspearman <- cor(kmat, method = "pearson")
diag(crspearman) <- NA

xlclipboard(crspearman)

table(tvc$status)

summary(abs(crspearman[lower.tri(crspearman)]))

fisher.test(table(kmat[,1], kmat[,2]))

cor(kmat[,1], kmat[,2])

stop()

#########################################################
# FULL MULTIVARIATE MODEL

if(rebuild) {
  cxme_full_id <- coxme(Surv(length, status) ~ 
                          itv_info_campaigns + itv_intl_travel + itv_school +
                          itv_public_events + itv_gatherings + itv_workplace + 
                          itv_movement + itv_public_transport + itv_stay_at_home + 
                          (1|id) + (1|country),
                        data = tvc)
  save(cxme_full_id, file = "rdata/cxme_full_id_200716.Rdata")
}

load("rdata/cxme_full_id_200716.Rdata")

coefmat <- coxme_coef( cxme_full_id )

coefmat

percs <- cbind(rownames(coefmat), to_percents(data.table(coefmat)))

# stop()

print(load("rdata/delay_table_200518.Rdata"))

svgf("fig/fig_3_delay_mv_coeffs", 12, 3.5)
{
  par(mfrow = c(1,2))
  # Panel, intervention delay
  {
    par(mar = c(5, 15, 2, 0))
    xlabtitle <- "Implementation delay from\n1st local divergence (days)"
    plot(0,0,type = "n", xlab = xlabtitle, ylab = "",
         xlim = c(-80, +80), ylim = c(0, length(itv_vars)) + 0.5, yaxt = "n", xaxt = "n")
    for(i in 1:length(itv_vars)) abline(h = i, lty = 2, col = "lightgrey")
    boxplot(delay ~ itv, delay_table, horizontal = TRUE, las = 1,
            xlab = xlabtitle, ylab = "", lwd = 2,
            col = rgb(0.6, 0, 0, 0.3), border = rgb(0.6, 0, 0), outline = FALSE, add = TRUE)
    beeswarm(delay ~ itv, delay_table, add = T, pch = 19, vertical = FALSE, col = rgb(0.6, 0, 0, 0.5), method = "hex", cex = 0.5)
  }
  # Panel, coefficients
  {
    par(mar = c(5, 1, 2, 20))
    plot(1, 0, type = "n", log = "", xlim = c(-55, 25), ylim = c(nrow(percs), 0) + 0.5, yaxt = "n",
         xlab = "Relative reproduction number\n(% change)", ylab = "")
    abline(v = 1, lty = 2, col = "lightgray")
    for(i in 1:nrow(percs)) {
      abline(h = i, lty = 2, col = "lightgrey")
    }
    points(percs$coef, 1:nrow(percs), col = rgb(0.6, 0, 0), lwd = 2, cex = 1.2)
    arrows(percs$lower, 1:nrow(percs), percs$upper, 1:nrow(percs),
           length = 0.02, angle = 90, code = 3, col = rgb(0.6, 0, 0), lwd = 2, lend = 2)
    # axis(2, at = 1:nrow(percs), labels = itv_labels, las = 2)
  }
  
}
dev.off()

#########################################################
# COVARIATE INFLUENCE MATRIX

# Drop covariates one at a time and measure influence

# Question of the use of nested vs crossed random effects
# https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified

if(rebuild) {
  library(parallel)
  
  cl <- makeCluster(length(itv_vars))
  
  clusterEvalQ(cl, {
    library(data.table)
    library(coxme)
  })
  
  clusterExport(cl, c("tvc", "cxme_full_id"))
  
  clusterApply(cl, itv_varlist, function(dropped) {
    save(dropped, file = sprintf("rdata/flag_dropped_%s_200716.Rdata", dropped))
    drop_model <- update(cxme_full_id, formula(
      sprintf("~ . - %s", dropped)
    ))
    save(drop_model, file = sprintf("rdata/cxme_dropped_%s_200716.Rdata", dropped))
  })
  
  stopCluster(cl)
}


model_list <- list()

for(dropped in itv_vars) {
  load(sprintf("rdata/cxme_dropped_%s_200716.Rdata", dropped))
  model_list[[dropped]] <- drop_model
}

for(dropped in itv_vars) {
  print(dropped)
  print(
    anova(cxme_full_id, model_list[[dropped]])
  )
}

model_coefs <- lapply(model_list, coxme_coef)

# Repeated merge
z <- data.table(var = rownames(coefmat), full = coefmat[,"coef"])
for(dropped in c(itv_vars)) {
  zz <- data.table(
    var = rownames(model_coefs[[dropped]]),
    coef = model_coefs[[dropped]][,"coef"])
  setnames(zz, "coef", dropped)
  z <- merge(z, zz, by = "var", all.x = TRUE, sort = FALSE)
}

z


zdiff <- data.table(z)
zdiff[ , (itv_vars) := lapply(.SD, function(x) x - full), .SDcols = itv_vars]

# Clean labels
zdiff[, var := itv_labels]
setnames(zdiff, itv_vars, itv_labels)

xlclipboard(zdiff)

zdiff

# Change to percents before taking diff
zperc <- data.table(z)
zperc[, full := 100*(exp(full) - 1)]
zperc[ , (itv_vars) := lapply(.SD, function(x) 100*(exp(x) - 1)), .SDcols = itv_vars]
zperc[ , (itv_vars) := lapply(.SD, function(x) x - full), .SDcols = itv_vars]
zperc[, var := itv_labels]
setnames(zperc, itv_vars, itv_labels)

zperc
xlclipboard(zperc)

#########################################################
# FULL MULTIVARIATE MODEL CONTROLLED ON TIME (TAXONOMIC BIAS)

if(rebuild) {
  cxme_full_id_time <- coxme(Surv(length, status) ~ time_source +
                               itv_info_campaigns + itv_intl_travel + itv_school +
                               itv_public_events + itv_gatherings + itv_workplace + 
                               itv_movement + itv_public_transport + itv_stay_at_home + 
                               (1|id) + (1|country),
                             data = tvc)
  save(cxme_full_id_time, file = "rdata/cxme_full_id_time_200716.Rdata")
}

load("rdata/cxme_full_id_time_200716.Rdata")

coefmat_time <- coxme_coef( cxme_full_id_time )

coefmat_time

# Express time influence in months
coefmat_time[1,1] <- coefmat_time[1,1] / 12
coefmat_time[1,3] <- coefmat_time[1,3] / 12

percs_time <- cbind(rownames(coefmat_time), to_percents(data.table(coefmat_time)))

#########################################################
# TABLE S3 - Model comparison with and without time

percs

percs_time

tables3 <- merge(percs_time, percs, by = "V1", all.x = TRUE, sort = FALSE)

tables3[, Factor := c("Elapsed time (per month)", itv_labels)]
tables3[, `Base model` := sprintf("%.1f (%.1f to %.1f)", coef.y, lower.y, upper.y)]
tables3[, `Time-adjusted model` := sprintf("%.1f (%.1f to %.1f)", coef.x, lower.x, upper.x)]

xlclipboard(tables3)


