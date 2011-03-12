plot.plspm <-
function(x, what="inner", how="joint", arr.pos=0.5, box.prop=.5, box.cex=1, cex.txt=1, ...)
{
    # ======================= plot function ============================
    # Displays path diagrams for each block of the outer model
    # =========================== ARGUMENTS ============================
    # x: object of class "plspm" obtained from function "plspm" 
    # what: a character string indicating what to plot
    #        options: "loadings", or "weights"
    # arr.pos: relative position of arrowhead on arrow segment/curve
    # box.prop: length/width ratio of label box
    # box.cex: relative size of text in boxes
    # cex.txt: relative size of arrow text

    # ==================== Checking function arguments =================
    WHAT <- c("inner","loadings","weights","all")
    what <- match.arg(what, WHAT)
    HOW <- c("joint","split")
    how <- match.arg(how, HOW)

    IDM <- x$model$IDM
    blocks <- x$model$blocks    
    modes <- x$model$modes
    lvs <- nrow(IDM)

    if (how=="joint")
    {
        rs <- c(1,1,1,2,2,2,2,2,3,2,3,3)
        cs <- c(1,2,3,2,3,3,4,4,3,5,4,4)
        index.mat <- cbind(1:12,rs,cs)
        if (what=="inner")
            .innerplot(x$path.coefs, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt)
        if (what=="loadings")
        {
            dev.new()
            par(mfrow=index.mat[lvs,2:3]) 
            par(mar=c(0,3,2.5,2)) 
            .loadingsplot(IDM, modes, blocks, x$loadings, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt, newdev=FALSE)
        }
        if (what=="weights")
        {
            dev.new()
            par(mfrow=index.mat[lvs,2:3]) 
            par(mar=c(0,3,2.5,2)) 
            .weightsplot(IDM, blocks, x$out.weights, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt, newdev=FALSE)
        }
        if (what=="all")
        {
            .innerplot(x$path.coefs, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt)
            dev.new()
            par(mfrow=index.mat[lvs,2:3]) 
            par(mar=c(0,3,2.5,2)) 
            .loadingsplot(IDM, modes, blocks, x$loadings, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt, newdev=FALSE)
            dev.new()
            par(mfrow=index.mat[lvs,2:3]) 
            par(mar=c(0,3,2.5,2)) 
            .weightsplot(IDM, blocks, x$out.weights, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt, newdev=FALSE)
        }
    } else { # how=="split"
        if (what=="inner")
            .innerplot(x$path.coefs, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt)
        if (what=="loadings")
            .loadingsplot(IDM, modes, blocks, x$loadings, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt, newdev=TRUE)
        if (what=="weights")
            .weightsplot(IDM, blocks, x$out.weights, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt, newdev=TRUE)
        if (what=="all")
        {
            .innerplot(x$path.coefs, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt)
            .loadingsplot(IDM, modes, blocks, x$loadings, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt, newdev=TRUE)
            .weightsplot(IDM, blocks, x$out.weights, arr.pos=arr.pos,
                 box.prop=box.prop, box.cex=box.cex, cex.txt=cex.txt, newdev=TRUE)
        }
    }

} # end function

