################################################################################
# 
# Project:  "enmpa: An R package for ecological niche modeling using
#            presence-absence data and generalized linear models"
# 
# Authors:   Luis F. Arias-Giraldo (lfarias.giraldo@gmail.com)
#            Marlon E. Cobos       (manubio13@gmail.com) 
# Date:      Feb 2024
#
# Purpose:   This script provides code to reproduce analyses and figures for the 
#            use of "enmpa" based on simulated data for a pathogen species
#
################################################################################

#~ R packages needed
library(enmpa)
library(terra)

################################################################################
#~ Load species occurrences and environmental data

data("enm_data", package = "enmpa")
vars <- rast(system.file("extdata", "vars.tif", package = "enmpa"))

head(enm_data)
plot(vars, nr = 2)

################################################################################
#~ Exploration of variables for models.
#~ Adapted methods developed by Cobos and Peterson (2022) to identify 
#~ relevant variables for characterizing species' ecological niches

#~ Univariate non-paramentric test
sn_uni_bio1  <- niche_signal(data = enm_data, variables = "bio_1",
                             condition = "Sp", method = "univariate")

sn_uni_bio12 <- niche_signal(data = enm_data, variables = "bio_12",
                             condition = "Sp", method = "univariate")

#~ Plotting results 
enmpa::plot_niche_signal(sn_uni_bio1, variables = "bio_1",
                         xlab = "BIO 1 (mean)")

enmpa::plot_niche_signal(sn_uni_bio12, variables = "bio_12",
                         xlab = "BIO 12 (mean)")


#~ Multivariate analysis based on a PERMANOVA
#~ task memory-intensive & time-consuming
sn_per  <- niche_signal(data = enm_data, variables = c("bio_1", "bio_12"),
                        condition = "Sp", method = "permanova")

sn_per

#~ plotting results
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)

enmpa::plot_niche_signal(sn_per, variables = c("bio_1", "bio_12"),
                         pch = c(16, 4),
                         xlab = "BIO 1", ylab = "BIO 12",
                         ellipses = FALSE)

legend("topright", legend = c("Host records", "Pathogen present"),
       cex = 0.9, pch = c(16, 4), pt.cex = c(1.5, 1.5), 
       col = c("black", "red"), 
       y.intersp = 2, bty = "n")

enmpa::plot_niche_signal(sn_per, variables = c("bio_1", "bio_12"),
                         xlab = "BIO 1", ylab = "BIO 12",
                         ellipses = TRUE)

pv <- sn_per$permanova_results$`Pr(>F)`[1]
pvf <- ifelse (pv < 0.05, "< 0.05", paste0( "= ",round(pv, 2)))

legend("topleft", legend = substitute(paste(italic("p "), x), env = list(x = pvf)),
       bty = "n", cex = 1, inset = 0.01)

legend("topright", legend = c("Host niche", "Pathogen niche"),  cex = 0.9,
       lty = 1, lwd = 1.5, col = c("black", "red"), 
       y.intersp = 2, bty = "n")

################################################################################
#~ Get formulas combination

#~ considering linear (l), quadratic (q) & products or two-way interactions (p)
get_formulas(dependent = "Sp", independent = c("bio_1", "bio_12"),
             type = "lpq", mode = "light")

get_formulas(dependent = "Sp", independent = c("bio_1", "bio_12"),
             type = "lpq", mode = "moderate")

get_formulas(dependent = "Sp", independent = c("bio_1", "bio_12"),
             type = "lpq", mode = "intensive")

################################################################################
#~ GLM model calibration
calibration <- calibration_glm(data = enm_data,
                               dependent = "Sp",
                               independent = c("bio_1", "bio_12"),
                               formula_mode = "intensive",
                               response_type = "lpq", 
                               out_dir = "calibration_output",
                               cv_kfolds = 5,
                               parallel = F,
                               exclude_bimodal = T)

calibration$calibration_results  # Results obtained from cross-validation for all models
calibration$summary              # A summary of statistics for all models
calibration$selected             # Selected models
calibration$data                 # Input data used

#~ Fitting selected models
fit_fmodes <- fit_selected(glm_calibration = calibration)
fit_fmodes

################################################################################
#~ Response curves

#~ Response curves for BIO-1 for models ID 29, 31 and consensus
par(mfrow = c(2, 3), mar = c(5, 4, 4, 2) + 0.1)

response_curve(fit_fmodes$ModelID_29, variable =  "bio_1",  xlab = "BIO-1",
               new_data = vars, main = "Model ID 29",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)

response_curve(fit_fmodes$ModelID_31, variable =  "bio_1",  xlab = "BIO-1",
               new_data = vars,  main = "Model ID 31",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)

response_curve(fit_fmodes, variable =  "bio_1",  xlab = "BIO-1",
               new_data = vars,  main = "Consensus",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)


#~ Response curves for BIO-12 for models ID 29, 31 and consensus

response_curve(fit_fmodes$ModelID_29, variable =  "bio_12",  xlab = "BIO-12",
               new_data = vars, main = "Model ID 29",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)

response_curve(fit_fmodes$ModelID_31, variable =  "bio_12",  xlab = "BIO-12",
               new_data = vars,  main = "Model ID 31",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)

response_curve(fit_fmodes, variable =  "bio_12",  xlab = "BIO-12",
               new_data = vars,  main = "Consensus",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)


################################################################################
#~ Variable importance based on explained deviance

(vi_29 <- var_importance(fit_fmodes$ModelID_29))
(vi_31 <- var_importance(fit_fmodes$ModelID_31))
(vi_c <- var_importance(fit_fmodes))

vi_29$predictor <- gsub("bio_", "BIO-", vi_29$predictor)
vi_31$predictor <- gsub("bio_", "BIO-", vi_31$predictor)
vi_c$predictor <- gsub("bio_", "BIO-", vi_c$predictor)

#~ Plots
dev.off()
enmpa::plot_importance(vi_29)
enmpa::plot_importance(vi_31)
enmpa::plot_importance(vi_c)

################################################################################
#~ Final model projection onto the area of North America

###~ Individual model predictions
preds <- predict_selected(fit_fmodes,
                          newdata = vars,
                          consensus = TRUE)

terra::plot(preds$predictions, nr = 2)
terra::plot(preds$consensus)

################################################################################
#~ Final model evaluation with an independent data

#~ Loading an independent dataset
data("test", package = "enmpa")
head(test)

#~ The independent evaluation data are divided into two groups: 
#~ presences-absences (test_01) and presences-only (test_1).
test_1 <- test[test$Sp == 1,]
test_01 <- test

#~ Final evaluation using the weighted average based on Akaike weights
wmean <- preds$consensus$Weighted_average

##~ Presences-absences evaluation
ie_01 <- independent_eval01(prediction = wmean,
                            observation = test_01$Sp,
                            lon_lat = test_01[, c("lon", "lat")])
ie_01


##~ Presences-only evaluation

#~ thresholds for only presences
t1 = ie_01[1,"Threshold"] # ESS
t2 = ie_01[2,"Threshold"] # maxTSS
t3 = ie_01[3,"Threshold"] # SEN90

independent_eval1(prediction = wmean,
                  threshold = t1,
                  lon_lat = test_1[, c("lon", "lat")])

independent_eval1(prediction = wmean,
                  threshold = t2,
                  lon_lat = test_1[, c("lon", "lat")])

independent_eval1(prediction = wmean,
                  threshold = t3,
                  lon_lat = test_1[, c("lon", "lat")])

#~ end of the script

