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
# BETWEEN-COUNTRY COMPARISON OF DIVERGENCE RATES
 
source("helpers.R")

print(load("rdata/nextstrain_tree_200512.Rdata"))
print(load("rdata/branchlengths_200512.Rdata"))
print(load("rdata/country_baseline_200512.Rdata"))
print(load("rdata/interventions_200518.Rdata"))
print(load("rdata/country_table_200518.Rdata"))

####################################################################
# Pool countries with <10 events into an "Others" category

# Minimal no. of internal branches (events) to be included in figure
min_events <- 10

tb <- table(blct[tip == FALSE]$country)
tb <- sort(tb, decreasing = T)

country_baseline[, country_with_others := factor(as.character(country), levels = names(tb))]
blct[ , country_with_others := factor(as.character(country), levels = names(tb))]
names(tb)[tb < min_events] <- "Others"
levels(blct$country_with_others) <- names(tb)
blct[ , country_with_others := relevel(country_with_others, "China")]
levels(blct$country_with_others)

####################################################################
# TIME SCALES Express length in days and delays in weeks

blct[, length := length * 364.25]
blct[, delay_global := delay_global / 7]
blct[, delay_in_country := delay_in_country / 7]

####################################################################
# FIGURE S2 - Correlation of first divergence dates with declared 10th case and 10th death

country_baseline <- merge(country_baseline, blct[ , .(first_div = first_div[1]), by = country], by = "country", all.x = T, sort = F)

Sys.setlocale("LC_TIME", "C") # Dates in C standard display
pdff("fig/cases_vs_divergence_dates", 10, 5)
{
  par(mfrow = c(1,2))
  par(mar = c(5,5,2,2))
  xl <- ymd("2019-12-01", "2020-05-15")
  ft <- country_baseline$country %in% levels(blct$country_with_others) & !is.na(country_baseline$date_10cases) & !is.na(country_baseline$date_10deaths)
  x_cases <- country_baseline[ ft ]$date_10cases
  x_deaths <- country_baseline[ ft ]$date_10deaths
  y <- dealigndate(country_baseline[ ft ]$first_div)
  labels <- country_baseline[ ft ]$country
  
  plot(x_cases, y, type = "n", xlim = xl, ylim = xl, pch = 19,
       xlab = "Date of 10th reported Covid-19 case", ylab = "Date of first detected divergence event")
  abline(0,1, lty = 2, col = "lightgray")
  textplot(x_cases, y, words = labels, new = FALSE, cex = 0.8)
  r2 <- cor(as.numeric(x_cases), as.numeric(y))
  legend("top", bty = "n", title = sprintf("Correlation = %.2f", r2), legend = "")
  
  plot(x_deaths, y, type = "n", xlim = xl, ylim = xl,
       xlab = "Date of 10th reported Covid-19 death", ylab = "Date of first detected divergence event")
  abline(0,1, lty = 2, col = "lightgray")
  textplot(x_deaths, y, words = labels, new = FALSE, cex = 0.8)
  r2 <- cor(as.numeric(x_deaths), as.numeric(y))
  legend("top", bty = "n", title = sprintf("Correlation = %.2f", r2), legend = "")
}
dev.off()



####################################################################
# FIGURE 2 - DIVERGENCE RATES PER COUNTRY

# Cox model
cx <- coxph(Surv(length, !tip) ~ country_with_others, blct, x = TRUE)
cx
cxz <- cox.zph(cx)
cxz
cxf <- coef(summary(cx))
cxf <- data.table(country = levels(blct$country_with_others)[-1], cxf, exp(confint(cx)))

setnames(cxf, c("country", "coef", "exp_coef", "se", "z", "p", "lo", "hi"))
cxf <- cxf[order(-coef)]
cxf

# Add China to align with others
cxf <- rbind(cxf, data.table(
  country = "China (reference)", coef = 0, exp_coef = 1, se = 0, z = 0, p = 0, lo = NA, hi = NA
))

cxf

###########################################################################
# FIGURE 2 - DISTRIBUTION OF DIVERGENCE EVENTS PER COUNTRY

# Add display name for DRC
cxf[country == "Democratic Republic of the Congo", country := "DR Congo"]

# Assign proper date to each branch source
blct[, display_date := dealigndate(time_source)]
blct[, display_country := country_with_others]
# Make a boxplot/beeswarm per country, aligned with the above figure
blct$display_country
levels(blct$display_country)[which(levels(blct$display_country) == "China")] <- "China (reference)"
levels(blct$display_country)[which(levels(blct$display_country) == "Democratic Republic of the Congo")] <- "DR Congo"
# Align with previous plot (trick here: can't reverse boxplot ordering using ylim, so reverse level ordering)
blct[, display_country := factor(as.character(display_country), levels = rev(cxf$country))]

# UPDATE 2020-07-26
# Use percent changes for consistency

cxf[, perc := 100*(exp_coef - 1)]
cxf[, perc_lo := 100*(lo - 1)]
cxf[, perc_hi := 100*(hi - 1)]


###########################################################################
# FIGURE 2 - GENERATE PANELS C AND D

svgf("fig/country_cox_time_c10_perc", 11, 7)
{
  par(mfrow = c(1, 2))
  n_countries <- nrow(cxf)
  
  # Plot divergence dates
  {
    par(mar = c(5, 10, 4, 0.5))
    
    par(lend = 2)
    plot(0, 0, xlim = range(blct$time_source), ylim = c(1, n_countries), xaxt = "n", yaxt = "n",
         xlab = "Inferred dates of transmission events\n", ylab = "")
    axis.at <- aligndate(ymd(c("2020-01-01", "2020-02-01", "2020-03-01", "2020-04-01")))
    axis(1, at = axis.at,
         labels = c("January", "February", "March", "April"))
    for(i in 1:nlevels(blct$display_country)) {
      abline(h = i, lty = 2, col = "lightgrey")
    }
    for(at in axis.at) {
      abline(v = at, lty = 2, col = "lightgrey")
    }
    boxplot(time_source ~ display_country, blct, horizontal = TRUE, las = 1,
            lwd = 2, xaxt = "n", yaxt = "n",
            col = rgb(0.6, 0, 0, 0.3), border = rgb(0.6, 0, 0), outline = TRUE, add = TRUE)
    axis(2, at = nrow(cxf):1, labels = cxf$country, las = 2)
  }
  
  # Plot divergence HR
  {
    par(mar = c(5, 0.5, 4, 10))

    plot(1, 0, type = "n", log = "", xlim = c(-80, +40), ylim = c(nrow(cxf), 1), yaxt = "n",
         xlab = "Relative reproduction number\n(% change)", ylab = "")
    abline(v = 0, lty = 2, col = "lightgray")
    for(i in 1:nlevels(blct$display_country)) {
      abline(h = i, lty = 2, col = "lightgrey")
    }
    points(cxf$perc, 1:nrow(cxf), col = rgb(0.6, 0, 0), lwd = 2, cex = 1.2)
    arrows(cxf$perc_lo, 1:nrow(cxf), cxf$perc_hi, 1:nrow(cxf),
           length = 0.02, angle = 90, code = 3, col = rgb(0.6, 0, 0), lwd = 2, lend = 2)
  }
}
dev.off()

############################################################
# FIGURE 2 - EXEMPLARY KAPLAN MEIER CURVES


svgf("fig/survplot_country", 8.5, 4)
{
  par(mfrow = c(1,2))
  par(mar = c(5,5,2,2))
  showct <- c("United States", "United Kingdom")
  svf <- survfit(Surv(length, !tip) ~ str,
                 blct[, .(length, tip, str = factor(as.character(country), levels = showct))])
  survplot(svf, "", showct)
  
  showct <- c("Spain", "Luxembourg")
  svf <- survfit(Surv(length, !tip) ~ str,
                 blct[, .(length, tip, str = factor(as.character(country), levels = showct))])
  survplot(svf, "", showct)
}
dev.off()

#############################################################
# FIGURE 2 - Correlation with case counts

names(country_baseline)
cxf

cxf2 <- merge(cxf, country_baseline, by = "country", all.x = TRUE, sort = FALSE)

# Dubious zero-cases/deaths etc
cxf2[ total_cases == 0, total_cases := NA ]
cxf2[ total_deaths == 0, total_deaths := NA ]
cxf2[ total_cases_per_million == 0, total_cases_per_million := NA ]
cxf2[ total_deaths_per_million == 0, total_deaths_per_million := NA ]

cor.test(cxf2$perc, log(cxf2$total_cases ), method = "pearson")
cor.test(cxf2$perc, log(cxf2$total_deaths ), method = "pearson")
cor.test(cxf2$perc, log(cxf2$total_cases_per_million ), method = "pearson")
cor.test(cxf2$perc, log(cxf2$total_deaths_per_million ), method = "pearson")

svgf("fig/casecounts_correlations", 6, 2.2)
{
  par(mfrow = c(1,2))
  par(mar = c(4, 6, 1, 2))
 
  options(scipen = 999) 
  plot(cxf2$perc, cxf2$total_cases, log = "y", las = 1, ylim = c(1, 2e6),
       xlab = "Relative reproduction number\n(% change relative to China)",
       ylab = "Counts (log scale)\n\n", pch = 19, col = rgb(0,0,0.8,0.7))
  points(cxf2$perc, cxf2$total_deaths, pch = 19, col = rgb(0.8,0,0,0.7))
  text(10, 100, "Cases", col = "blue", adj = 1)
  text(10, 10, "Deaths", col = "Red", adj = 1)
  
  plot(cxf2$perc, cxf2$total_cases_per_million, log = "y", las = 1, ylim = c(1e-1, 1e4),
       xlab = "Relative reproduction number\n(% change relative to China)",
       ylab = "Counts per million\n(log scale)\n", pch = 19, col = rgb(0,0,0.8,0.7))
  points(cxf2$perc, cxf2$total_deaths_per_million, pch = 19, col = rgb(0.8,0,0,0.7))
  text(10, 10, "Cases", col = "blue", adj = 1)
  text(10, 1, "Deaths", col = "Red", adj = 1)
}
dev.off()


plot(cxf2$perc, cxf2$total_cases_per_million, log = "y")
points(cxf2$perc, cxf2$total_deaths_per_million)
