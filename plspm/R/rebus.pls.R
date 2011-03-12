rebus.pls <-
function(pls, Y=NULL, stop.crit=0.005, iter.max=100)
{
    # ========================= rebus.pls function ========================
    # Function to perform Response-Based Unit Segmentation (REBUS)
    # in a single function. It integrates "res.clus" and "rebus"
    # =========================== ARGUMENTS ===============================
    # pls: object of class "plspm"
    # Y: optional data matrix used when pls$data is null
    # stop.crit: stop criterion number (must be between 0 and 1)
    # iter.max: maximum number of iterations (must be an integer)

    # ==================== Checking function arguments ===================
    if (class(pls)!="plspm") 
        stop("argument 'pls' must be an object of class 'plspm'")
    if (any(pls$model$modes!="A"))# checking reflective modes
        stop("REBUS only works for reflective modes")
    if (!pls$model$scaled)# checking scaled data
        stop("REBUS only works with scaled='TRUE'")
    if (!is.null(Y)) # if Y available
    {
        if (is.null(pls$data))
        {
            if (!is.matrix(Y) && !is.data.frame(Y))
                stop("Invalid object 'Y'. Must be a numeric matrix or data frame.")
            if (nrow(Y)!=nrow(pls$latents))
                stop("Argument 'pls' and 'Y' are incompatible. Different number of rows.")
        }
    } else { # if no Y
        if (is.null(pls$data)) 
            stop("Argument 'Y' is missing. No dataset available.")
    }
    if (mode(stop.crit)!="numeric" || length(stop.crit)!=1 || 
        stop.crit<0 || stop.crit>=1)
    {
        warning("Invalid stop criterion 'stop.crit'. Deafult value 0.005 is used")
        stop.crit <- 0.005
    }
    if (mode(iter.max)!="numeric" || length(iter.max)!=1 || 
        iter.max<=1 || (iter.max%%1)!=0)
    {
        warning("Invalid number of maximum iterations 'iter.max'. Deafult value 100 is used")
        iter.max <- 100
    }

    resid <- res.clus(pls, Y)
    print("Please enter the number of classes (an integer > 1), and then press Enter:")
    nk <- scan(what=integer(1), n=1)
    if (mode(nk)!="numeric" || length(nk)!=1 || 
        nk<=1 || (nk%%1)!=0)
        stop("Invalid number of classes. Must be an integer larger than 1")    
    resul <- it.reb(pls, resid, nk, Y, stop.crit, iter.max)
    return(resul)
}

