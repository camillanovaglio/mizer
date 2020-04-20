##' @title Self-correlated serie
##' @param years Years of autocorrelated values
##' @param sigmaz error variance of the serie
##' @param phi autocorrelation lag coefficient
##' @param seed Initial Value for the seed in R. Just for calculation
##' @param fix.init Option to use an initial value
##' @return Residual autocorrelation or autocorrelation values
##' @author Javier Porobic
aut.rd <- function(years, sigmaz, phi, seed = NULL, fix.init = NULL){
    if(!is.null(seed)) set.seed(seed)
    rd     <- rnorm(years, sd = 3)
    sigmae <- sigmaz * (1 - phi ^ 2)
    zt     <- vector(mode = 'numeric',  length = years)
    zt[1]  <- sqrt(sigmaz) * rd[1]
    if(!is.null(fix.init)) zt[1] <- fix.init
    for(dt in 2 : years){
        zt[dt] <- zt[dt - 1] * phi + rd[dt] * sqrt(sigmae)
    }
    return(zt)
}




plot(1 : 100, aut.rd(years = 100, sigmaz = 0.9, phi = 0.99, fix.init = 10), type = 'l')
for(i in 1 : 4){
    lines(1 : 100, aut.rd(years = 100, sigmaz = 0.9, phi = 0.99, fix.init = 10), col = i + 1)
}
