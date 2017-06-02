coreNum = 14 # change this for the number of cores you have access to

lambdaRange <- seq(from = 0.1, to = 1, by = 0.1)
chronosModels <- c('correlated', 'relaxed', 'discrete')
chronosTest <- list(correlated = mclapply(lambdaRange, function(x) chronos(tr.uniques[[1]]$sppOnly,
                                          #calibration = tr.uniques[[1]]$sppOnly.calib,
                                          model = chronosModels[1], lambda = x), mc.cores = coreNum),
                    relaxed = mclapply(lambdaRange, function(x) chronos(tr.uniques[[1]]$sppOnly,
                                          #calibration = tr.uniques[[1]]$sppOnly.calib,
                                          model = chronosModels[2], lambda = x), mc.cores = coreNum),
                    discrete = mclapply(lambdaRange, function(x) chronos(tr.uniques[[1]]$sppOnly,
                                          #calibration = tr.uniques[[1]]$sppOnly.calib,
                                          model = chronosModels[3], lambda = x), mc.cores = coreNum))

chronosTest.mat <- matrix(NA, 10, 3, dimnames = list(as.character(lambdaRange), chronosModels))
for(i in seq(length(lambdaRange))) {
  for(j in chronosModels) {
    chronosTest.mat[i, j] <- attr(chronosTest[[j]][[i]], 'PHIIC')$PHIIC
  }
}
