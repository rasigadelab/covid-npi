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
# TIMING OF NON-PHARMACEUTICAL INTERVENTIONS

source("helpers.R")

print(load("rdata/nextstrain_tree_200512.Rdata"))
print(load("rdata/branchlengths_200512.Rdata"))
print(load("rdata/country_baseline_200512.Rdata"))
print(load("rdata/interventions_200518.Rdata"))
print(load("rdata/country_table_200518.Rdata"))

###########################################################################################
# DELAY FROM FIRST DIVERGENCE TO NPI

delay_table <- data.table(country_table)[, -c("n_genomes", "max_stringency")]
names(delay_table)

delay_table[, (itv_vars) := lapply(.SD, function(x) as.numeric(x - first_edge)), .SDcols = itv_vars][, -c("first_edge")]

# Melt
delay_table <- melt(delay_table, id.vars = "country", measure.vars = itv_vars, variable.name = "itv", value.name = "delay")
delay_table <- delay_table[is.finite(delay)]

summary(delay_table$delay)

# # Relabel
cbind(levels(delay_table$itv), itv_labels)
levels(delay_table$itv) <- itv_labels

# Reorder
delay_table[, itv := factor(as.character(itv), levels = c(
  rev(levels(itv))
))]

# Timing data will be used for SIR models
save(delay_table, file = "rdata/delay_table_200518.Rdata")

###########################################################################################
# FIGURE 3 - DELAY FROM FIRST DIVERGENCE TO NPI

# stop()

svgf("fig/intervention_delay", 6, 3.5)
{
  par(mfrow = c(1,1))
  par(mar = c(5, 15, 2, 2))
  # boxplot(delay ~ itv, delay_table, horizontal = TRUE, las = 1,
  #         xlab = "Delay from 1st divergence event to intervention (days)", ylab = "", lwd = 2,
  #         col = rgb(0.6, 0, 0, 0.3), border = rgb(0.6, 0, 0), outline = FALSE)
  plot(0,0,type = "n", xlab = "Delay from 1st divergence event\nto intervention (days)", ylab = "",
       xlim = c(-80, +80), ylim = c(0, length(itv_vars)) + 0.5, yaxt = "n", xaxt = "n")
  for(i in 1:length(itv_vars)) abline(h = i, lty = 2, col = "lightgrey")
  boxplot(delay ~ itv, delay_table, horizontal = TRUE, las = 1,
          xlab = "Delay from 1st divergence event to intervention (days)", ylab = "", lwd = 2,
          col = rgb(0.6, 0, 0, 0.3), border = rgb(0.6, 0, 0), outline = FALSE, add = TRUE)
  beeswarm(delay ~ itv, delay_table, add = T, pch = 19, vertical = FALSE, col = rgb(0.6, 0, 0, 0.5), method = "hex", cex = 0.5)
  
}
dev.off()

# Median delays
delay_table[, .(median = median(delay), mean = mean(delay), q1 = quantile(delay, 0.25), q3 = quantile(delay, 0.75)), by = itv]

##############################################################
# FIGURE S3 A - CONDITIONAL NPI IMPLEMENTATIONS

itv_between <- data.table(country_table)[, -c("n_genomes", "max_stringency")]
# Melt
itv_between <- melt(itv_between, id.vars = "country", measure.vars = itv_vars, variable.name = "itv", value.name = "date", na.rm = FALSE)
cbind(levels(itv_between$itv), itv_labels)
levels(itv_between$itv) <- itv_labels

# Cartesian product of intervention pairs by country
itv_between <- merge(itv_between, itv_between, by = "country", all.x = TRUE, all.y = TRUE, allow.cartesian = TRUE)
# Delay between interventions in each country (contains Inf when one date is missing, NaN when both are missing)
itv_between[, delay := date.y - date.x]
# Absolute delay
itv_between[, delay_abs := as.numeric((delay))][is.nan(delay), delay_abs := NA]

# Proportion of countries implementing each intervention
prop_itv <- itv_between[ itv.x == itv.y, .(n = sum(!is.na(delay_abs))), by = itv.x]

# Pretty print proportions for table/figure
n_ct <- 57
prop_itv[ , pret := sprintf("%g (%.1f)", n, 100 * n/n_ct)]

xlclipboard(prop_itv)

# SUmmarise absolute delays
itv_between_summary <- itv_between[ , .(
  q0 = min(delay_abs, na.rm = TRUE),
  q1 = quantile(delay_abs, 0.25, na.rm = TRUE),
  q2 = quantile(delay_abs, 0.50, na.rm = TRUE),
  q3 = quantile(delay_abs, 0.75, na.rm = TRUE),
  q4 = max(delay_abs, na.rm = TRUE),
  mean = mean(delay_abs, na.rm = TRUE),
  sd = sd(delay_abs, na.rm = TRUE),
  y_if_x = mean(delay < Inf, na.rm = T),
  both = mean(is.finite(delay_abs))
), by = .(itv.x, itv.y)]

itv_between_summary

# Matrix of conditional implementations
itv_between_mat <- dcast(itv_between_summary, itv.x ~ itv.y, fill = NA, value.var = "y_if_x")
itv_between_mat <- as.matrix(itv_between_mat[, -c("itv.x")])
rownames(itv_between_mat) <- colnames(itv_between_mat)
itv_between_mat

xlclipboard(itv_between_mat)

##############################################################
# FIGURE S3 B -  between-intervention delay (when both were implemented)
itv_between_summary_finite <- itv_between[ is.finite(delay_abs), .(
  q0 = min(delay_abs),
  q1 = quantile(delay_abs, 0.25),
  q2 = quantile(delay_abs, 0.50),
  q3 = quantile(delay_abs, 0.75),
  q4 = max(delay_abs),
  mean = mean(delay_abs),
  sd = sd(delay_abs),
  y_if_x = mean(delay < Inf, na.rm = T)
), by = .(itv.x, itv.y)]

# Delay from intervention in row to intervention in column when both interventions were implemented
itv_delay_mat <- dcast(itv_between_summary_finite, itv.x ~ itv.y, fill = NA, value.var = "q2")
itv_delay_mat <- as.matrix(itv_delay_mat[, -c("itv.x")])
rownames(itv_delay_mat) <- colnames(itv_delay_mat)
itv_delay_mat

xlclipboard(itv_delay_mat)
