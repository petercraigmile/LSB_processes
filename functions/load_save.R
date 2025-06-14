
save.it <- function (obj) {

    save(list=obj, file=paste("results/", obj, ".RData", sep=""))
}

load.it <- function (obj) {

    load(file=paste("results/", obj, ".RData", sep=""), envir=.GlobalEnv)
}
