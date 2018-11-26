
fitMx <- stan(file="Modelx.stan", data = data_multi, iter = 5000, chains = 4, control = list(adapt_delta = 0.95, stepsize_jitter = 0.5, max_treedepth = 13))

fitMx.dso <- new("cxxdso")
saveRDS(fitMx, file = "fitMx.rds")

