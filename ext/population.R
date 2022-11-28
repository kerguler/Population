dyn.load("./model.dylib")

numpar <- as.integer(0)
nummet <- as.integer(0)
tmp <- .C("init", numpar, nummet)
numpar <- tmp[[1]]
nummet <- tmp[[2]]



result <- rep(0, finalT * ArboRisk::parameters$internal$vector_columns$NumMet)

