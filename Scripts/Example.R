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
#~ to get the last dev version: "remotes::install_github("Luisagi/enmpa")" 
library(enmpa) 
library(terra)

#~ set terrain default palette
options(terra.pal=terrain.colors(234, rev = TRUE))

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

#~ Multivariate analysis based on a PERMANOVA___________________________________
#~ task memory-intensive & time-consuming
sn_per  <- niche_signal(data = enm_data, variables = c("bio_1", "bio_12"),
                        condition = "Sp", method = "permanova")

sn_per

#~ plotting results: Figure S1

#png(filename = "Figure S1.png", width = 300, height = 180, units = "mm", res = 300)
par(mfrow = c(1, 2), mar = c(5, 4, 2, 1) + 0.1)

enmpa::plot_niche_signal(sn_per, variables = c("bio_1", "bio_12"),
                         pch = c(16, 4),
                         xlab = "BIO 1", ylab = "BIO 12",
                         ellipses = FALSE)

legend("topright", legend = c("Host records", "Pathogen present"),
       cex = 0.9, pch = c(16, 4), pt.cex = c(1, 1), 
       col = c("black", "red"), 
       y.intersp = 1, bty = "n", inset = c(0,0))

enmpa::plot_niche_signal(sn_per, variables = c("bio_1", "bio_12"),
                         xlab = "BIO 1", ylab = "BIO 12",
                         ellipses = TRUE)

pv <- sn_per$permanova_results$`Pr(>F)`[1]
pvf <- ifelse (pv < 0.05, "< 0.05", paste0( "= ", round(pv, 2)))

legend("topleft", legend = substitute(paste(italic("p "), x), env = list(x = pvf)),
       bty = "n", cex = 1, inset = c(0,0))

legend("topright", legend = c("Host niche", "Pathogen niche"),  cex = 0.8,
       lty = 1, lwd = 1.5, col = c("black", "red"), 
       y.intersp = 1, bty = "n", inset = c(0,0))

#dev.off()


#~ Univariate non-parametric test _____________________________________________
sn_uni_bio1  <- niche_signal(data = enm_data, variables = "bio_1",
                             condition = "Sp", method = "univariate")

sn_uni_bio12 <- niche_signal(data = enm_data, variables = "bio_12",
                             condition = "Sp", method = "univariate")

#~ Plotting results: Figure S2

#png(filename = "Figure S2.png", width = 300, height = 180, units = "mm", res = 300)
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1)
enmpa::plot_niche_signal(sn_uni_bio1, variables = "bio_1",  statistic = "mean",
                         xlab = "BIO 1 (mean)")

enmpa::plot_niche_signal(sn_uni_bio12, variables = "bio_12", statistic = "mean",
                         xlab = "BIO 12 (mean)")

enmpa::plot_niche_signal(sn_uni_bio1, variables = "bio_1",  statistic = "SD",
                         xlab = "BIO 1 (SD)")

enmpa::plot_niche_signal(sn_uni_bio12, variables = "bio_12", statistic = "SD",
                         xlab = "BIO 12 (SD)")
#dev.off()

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
                               exclude_bimodal = T, 
                               seed = 1)
calibration
summary(calibration)

calibration$calibration_results  # Results obtained from cross-validation for all models
calibration$summary              # A summary of statistics for all models
calibration$selected             # Selected models
calibration$data                 # Input data used

#~ Fitting selected models
fmodels <- fit_selected(glm_calibration = calibration)
summary(fmodels)


################################################################################
#~ Response curves (Figure 1)

#~ Response curves for BIO-1 for models ID 29, 31 and consensus

#png(filename = "Figure 1A.png", width = 300, height = 180, units = "mm", res = 300)
par(mfrow = c(2, 3), mar = c(5, 4, 4, 2) + 0.1)

response_curve(fitted = fmodels, modelID = "ModelID_29",
               variable =  "bio_1",  xlab = "BIO-1",
               new_data = vars, main = "Model ID 29",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)

response_curve(fitted = fmodels, modelID = "ModelID_31",
               variable =  "bio_1",  xlab = "BIO-1",
               new_data = vars,  main = "Model ID 31",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)

response_curve(fmodels, variable =  "bio_1",  xlab = "BIO-1",
               new_data = vars,  main = "Consensus",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)


#~ Response curves for BIO-12 for models ID 29, 31 and consensus

response_curve(fitted = fmodels, modelID = "ModelID_29", 
               variable =  "bio_12",  xlab = "BIO-12",
               new_data = vars, main = "Model ID 29",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)

response_curve(fitted = fmodels, modelID = "ModelID_31",
               variable =  "bio_12",  xlab = "BIO-12",
               new_data = vars,  main = "Model ID 31",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)

response_curve(fmodels, variable =  "bio_12",  xlab = "BIO-12",
               new_data = vars,  main = "Consensus",
               cex.lab=1.25, cex.axis=1.15, cex.main=1.25, cex.sub=1.25)

dev.off()

################################################################################
#~ Variable importance based on explained deviance

(vi_29 <- var_importance(fitted = fmodels, modelID = "ModelID_29"))
(vi_31 <- var_importance(fitted = fmodels, modelID = "ModelID_31"))
(vi_c <- var_importance(fmodels))

vi_29$predictor <- gsub("bio_", "BIO-", vi_29$predictor)
vi_31$predictor <- gsub("bio_", "BIO-", vi_31$predictor)
vi_c$predictor <- gsub("bio_", "BIO-", vi_c$predictor)

#~ Plots
dev.off()
enmpa::plot_importance(vi_29)
enmpa::plot_importance(vi_31)

#png(filename = "Figure 1B.png", width = 300, height = 180, units = "mm", res = 300)
enmpa::plot_importance(vi_c, main = "Variable importance")
#dev.off()


################################################################################
#~ Final model projection onto the area of North America

###~ Individual model predictions
preds_E  <- predict_selected(fmodels, newdata = vars, consensus = TRUE,
                             extrapolation_type = "E")
preds_NE <- predict_selected(fmodels, newdata = vars, consensus = TRUE, 
                             extrapolation_type = "NE")
preds_EC <- predict_selected(fmodels, newdata = vars, consensus = TRUE,
                             extrapolation_type = "EC")
preds_E
preds_NE
preds_EC

terra::plot(c(preds_E$predictions,
              preds_NE$predictions,
              preds_EC$predictions),
            main = c("Model ID 29 (E)", "Model ID 31 (E)",
                     "Model ID 29 (NE)", "Model ID 31 (NE)",
                     "Model ID 29 (EC)", "Model ID 31 (EC)"),
            mar = c(1, 1, 2, 5), nr = 3)


###~ Consensus model predictions
terra::plot(preds_E$consensus, mar = c(4, 2, 2, 3))
terra::plot(preds_NE$consensus, mar = c(4, 2, 2, 3))
terra::plot(preds_EC$consensus, mar = c(4, 2, 2, 3))


#~ Figure S3
plot(c(preds_NE$predictions$ModelID_29,
       preds_EC$predictions$ModelID_29,
       preds_NE$predictions$ModelID_31,
       preds_EC$predictions$ModelID_31,
       preds_NE$consensus$Weighted_average,
       preds_EC$consensus$Weighted_average),
     main = c("Model ID 29 (NE)","Model ID 29 (EC)",
              "Model ID 31 (NE)", "Model ID 31 (EC)", 
              "Weigthed average (NE)", "Weigthed average (EC)"),
     mar = c(2.5, 2, 2, 5),
     nr = 3)

################################################################################
#~ Final model evaluation with an independent data

#~ Loading an independent dataset
data("test", package = "enmpa")
head(test)

#~ The independent evaluation data are divided into two groups: 
#~ presences-absences (test_01) and presences-only (test_1).
test_1 <- test[test$Sp == 1,]
test_01 <- test

##~ Presences-absences evaluation

projections <- list(
  ModelID_29 = preds$predictions$ModelID_29,
  ModelID_31 = preds$predictions$ModelID_31,
  Consensus_WA = preds$consensus$Weighted_average
)

ie_01 <- lapply(projects, function(x){
  independent_eval01(prediction = x,
                     observation = test_01$Sp,
                     lon_lat = test_01[, c("lon", "lat")])
})

ie_01

##~ Presences-only evaluation

#~ thresholds for each model using only presences
 
# ModelID_29
ModelID_29_T1 <- ie_01[["ModelID_29"]][1, "Threshold"] # ESS
ModelID_29_T2 <- ie_01[["ModelID_29"]][2, "Threshold"] # maxTSS
ModelID_29_T3 <- ie_01[["ModelID_29"]][3, "Threshold"] # SEN90

aux_list1 <- list(ESS = ModelID_29_T1, maxTSS = ModelID_29_T2, SEN90 = ModelID_29_T3)

lapply(aux_list1, function(th){
  independent_eval1(prediction = projects[["ModelID_29"]],
                    threshold = th,
                    lon_lat = test_1[, c("lon", "lat")])
  
})

# ModelID_31
ModelID_31_T1 <- ie_01[["ModelID_31"]][1, "Threshold"] # ESS
ModelID_31_T2 <- ie_01[["ModelID_31"]][2, "Threshold"] # maxTSS
ModelID_31_T3 <- ie_01[["ModelID_31"]][3, "Threshold"] # SEN90

aux_list2 <- list(ESS = ModelID_31_T1, maxTSS = ModelID_31_T2, SEN90 = ModelID_31_T3)

lapply(aux_list2, function(th){
  independent_eval1(prediction = projects[["ModelID_31"]],
                    threshold = th,
                    lon_lat = test_1[, c("lon", "lat")])
  
})


# Consensus_WA
Consensus_WA_T1 <- ie_01[["Consensus_WA"]][1, "Threshold"] # ESS
Consensus_WA_T2 <- ie_01[["Consensus_WA"]][2, "Threshold"] # maxTSS
Consensus_WA_T3 <- ie_01[["Consensus_WA"]][3, "Threshold"] # SEN90

aux_list3 <- list(ESS = Consensus_WA_T1, maxTSS = Consensus_WA_T2, SEN90 = Consensus_WA_T3)

lapply(aux_list3, function(th){
  independent_eval1(prediction = projects[["Consensus_WA"]],
                    threshold = th,
                    lon_lat = test_1[, c("lon", "lat")])
  
})

#~ end of the script

