
#################################### CALIBRAION ####################################

data_extraction_multiexperiment("Calibration_4")
fitCal4 <- stan(file="Model3.stan", data = data_multi, iter = 3000, warmup = 1000, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitCal4.dso <- new("cxxdso")
saveRDS(fitCal4, file = "fitCal4.rds")

data_extraction_multiexperiment("Calibration_5")
fitCal5 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitCal5.dso <- new("cxxdso")
saveRDS(fitCal5, file = "fitCal5.rds")

data_extraction_multiexperiment("Calibration_6")
fitCal6 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitCal6.dso <- new("cxxdso")
saveRDS(fitCal6, file = "fitCal6.rds")

#################################### DYNSTIM ####################################

data_extraction_multiexperiment("DynStim_1")
fitDyn1 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn1.dso <- new("cxxdso")
saveRDS(fitDyn1, file = "fitDyn1.rds")

data_extraction_multiexperiment("DynStim_2")
fitDyn2 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn2.dso <- new("cxxdso")
saveRDS(fitDyn2, file = "fitDyn2.rds")

data_extraction_multiexperiment("DynStim_3")
fitDyn3 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn3.dso <- new("cxxdso")
saveRDS(fitDyn3, file = "fitDyn3.rds")

data_extraction_multiexperiment("DynStim_4")
fitDyn4 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn4.dso <- new("cxxdso")
saveRDS(fitDyn4, file = "fitDyn4.rds")

data_extraction_multiexperiment("DynStim_5")
fitDyn5 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn5.dso <- new("cxxdso")
saveRDS(fitDyn5, file = "fitDyn5.rds")

data_extraction_multiexperiment("DynStim_6")
fitDyn6 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn6.dso <- new("cxxdso")
saveRDS(fitDyn6, file = "fitDyn6.rds")

data_extraction_multiexperiment("DynStim_7")
fitDyn7 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn7.dso <- new("cxxdso")
saveRDS(fitDyn7, file = "fitDyn7.rds")

data_extraction_multiexperiment("DynStim_8")
fitDyn8 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn8.dso <- new("cxxdso")
saveRDS(fitDyn8, file = "fitDyn8.rds")

data_extraction_multiexperiment("DynStim_9")
fitDyn9 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn9.dso <- new("cxxdso")
saveRDS(fitDyn9, file = "fitDyn9.rds")

data_extraction_multiexperiment("DynStim_10")
fitDyn10 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn10.dso <- new("cxxdso")
saveRDS(fitDyn10, file = "fitDyn10.rds")

data_extraction_multiexperiment("DynStim_11")
fitDyn11 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn11.dso <- new("cxxdso")
saveRDS(fitDyn11, file = "fitDyn11.rds")

data_extraction_multiexperiment("DynStim_12")
fitDyn12 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn12.dso <- new("cxxdso")
saveRDS(fitDyn12, file = "fitDyn12.rds")

data_extraction_multiexperiment("DynStim_13")
fitDyn13 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn13.dso <- new("cxxdso")
saveRDS(fitDyn13, file = "fitDyn13.rds")

data_extraction_multiexperiment("DynStim_14")
fitDyn14 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDyn14.dso <- new("cxxdso")
saveRDS(fitDyn14, file = "fitDyn14.rds")

data_extraction_multiexperiment("DynStimTooFast")
fitDynTooFast <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitDynTooFast.dso <- new("cxxdso")
saveRDS(fitDynTooFast, file = "fitDynTooFast.rds")

#################################### BANGBANG ####################################

data_extraction_multiexperiment("BangBang_1")
fitBang1 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitBang1.dso <- new("cxxdso")
saveRDS(fitBang1, file = "fitBang1.rds")

data_extraction_multiexperiment("BangBang_2")
fitBang2 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitBang2.dso <- new("cxxdso")
saveRDS(fitBang2, file = "fitBang2.rds")


#################################### PI ####################################

data_extraction_multiexperiment("PI_1")
fitPI1 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitPI1.dso <- new("cxxdso")
saveRDS(fitPI1, file = "fitPI1.rds")

data_extraction_multiexperiment("PI_2")
fitPI2 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitPI2.dso <- new("cxxdso")
saveRDS(fitPI2, file = "fitPI2.rds")

data_extraction_multiexperiment("PI_3")
fitPI3 <- stan(file="Model3.stan", data = data_multi, iter = 500, warmup = 450, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))
fitPI3.dso <- new("cxxdso")
saveRDS(fitPI3, file = "fitPI3.rds")










