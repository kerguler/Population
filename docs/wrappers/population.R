model <- function(filename) {
    ret <- list()

    dyn.load(filename)

    numpar <- as.integer(0)
    nummet <- as.integer(0)
    tmp <- .C("init", numpar, nummet)
    ret$numpar <- tmp[[1]]
    ret$nummet <- tmp[[2]]

    mpnames <- rep("", ret$numpar + ret$nummet)
    tmp <- .C("parnames", mpnames)
    ret$metnames <- tmp[[1]][1:ret$nummet]
    ret$parnames <- tmp[[1]][ret$nummet + 1:ret$numpar]

    ret$metids <- list(); for (i in 1:length(ret$metnames)) { ret$metids[ret$metnames[i]] <- i; }
    ret$parids <- list(); for (i in 1:length(ret$parnames)) { ret$parids[ret$parnames[i]] <- i; }

    sim <- function(envir, pr, ftime, rept=1) {
        envir <- as.double(envir)
        pr <- as.double(pr)
        ftime <- as.integer(ftime)
        rept <- as.integer(rept)
        res <- rep(0.0, rept * ftime * ret$nummet)
        success <- as.integer(0)
        tmp <- .C("sim",
                  envir,
                  pr,
                  ftime,
                  rept,
                  res,
                  success)
        res <- array(tmp[[5]], dim=c(ret$nummet, ftime, rept))
        return(res)
    }

    ret$sim <- sim

    return(ret)
}
