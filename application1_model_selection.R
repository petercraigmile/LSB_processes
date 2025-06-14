
## ======================================================================
## Script: application1_model_selection.R
## Purpose: Carry out model selection using NIC and BIC criterion for
## EEG series
## ======================================================================

library(spectral)
library(fields)
library(mvnfast)

source("functions/basis_functions.R")
source("functions/LSB_ARp.R")
source("functions/LSB_AR_SDF.R")

eeg <- scan("EEG_data/eeg_F3.txt")

T. <- length(eeg)

ts <- seq(from=0, len=length(eeg), by=1/256)

fs <- seq(0, 1/2, length=256)

p_values = 8:20
b_values = 2:8

param_grid = expand.grid(ar.basis.order = b_values,
                         ar.order = p_values)

p_grid = param_grid$ar.order
b_grid = param_grid$ar.basis.order

RNGkind("L'Ecuyer-CMRG")
set.seed(1)
parallel::mc.reset.stream()


# Run LSB-AR(p) models with X basis order = b and sigma.sq basis order = b
models = pbmcapply::pbmcmapply(function(p, b) 
  LSB.ARp.model.sel(eeg, p, b),
  p_grid,
  b_grid,
  mc.cores = parallel::detectCores() - 1)


nic.save = matrix(sapply(1:ncol(models), function(i) models["NIC",][[i]]),
                  nrow = length(p_values), ncol = length(b_values),
                  byrow = TRUE)
bic.save = matrix(sapply(1:ncol(models), function(i) models["BIC",][[i]]),
                  nrow = length(p_values), ncol = length(b_values),
                  byrow = TRUE)

## Selected model via NIC and BIC

nic.model <- models[c("ar.order", "basis.order"), which.min(models["NIC", ])]
bic.model <- models[c("ar.order", "basis.order"), which.min(models["BIC", ])]

# Choose order of basis for the time varying variance curve based on selected 
# model using NIC
nic.var.models <- pbmcapply::pbmclapply(b_values, function(s)
  LSB.ARp.model.sel(eeg, nic.model$ar.order,
                    nic.model$basis.order, s), 
  mc.cores = parallel::detectCores() - 1)

nic.var.basis.order <-
  nic.var.models[[which.min(sapply(1:length(nic.var.models), function(x)
    nic.var.models[[x]]$NIC))]]$var.basis.order

# Choose order of basis for the time varying variance curve based on selected 
# model using BIC

bic.var.models <- pbmcapply::pbmclapply(b_values, function(s)
  LSB.ARp.model.sel(eeg, bic.model$ar.order,
                    bic.model$basis.order, s), 
  mc.cores = parallel::detectCores() - 1)

bic.var.basis.order <-
  bic.var.models[[which.min(sapply(1:length(bic.var.models), function(x)
    bic.var.models[[x]]$BIC))]]$var.basis.order

save(models, nic.var.models, bic.var.models,
     file = "results/model_selection_full_updated.RData")

save(nic.model, bic.model, nic.var.basis.order, bic.var.basis.order,
     file = paste("results/models_orders.RData"))
