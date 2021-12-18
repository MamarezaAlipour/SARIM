SarimStartupMessage <- function()
{
    # Startup message obtained as
    # > figlet -f slant MCLUST
    msg <- c(paste0(
        packageVersion("Sarim")),
        "\n  Structured Additive Regression Models using Iterative Methods.")
    return(msg)
}

.onAttach <- function(lib, pkg)
{
    # unlock .mclust variable allowing its modification
    ## unlockBinding("sarim", asNamespace("sarim"))
    # startup message
    msg <- SarimStartupMessage()
    if (!interactive())
        msg[1] <- paste("Package 'Sarim' version", packageVersion("Sarim"))
    packageStartupMessage(msg)
    invisible()
}