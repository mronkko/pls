resclus.plot <-
function(x, k=2, col.up="black", col.down=rainbow(k), lty.up=2, lty.down=1,
                    lwd.up=1, lwd.down=1, type="rectangle", knot.pos="mean",
                    show.labels=FALSE, only.tree=FALSE, members, ...)
{
    # ======================= resclus.plot function ========================
    # Displays a colored dendrogram for REBUS partitions (res.clus function)
    # =========================== ARGUMENTS ================================
    # x: object of class "hclust" obtained from function "res.clus" 
    # k: the number of classes
    # col.up: The color to be used for the lines above the cut of the dendrogram
    # col.down: The colors to be used for the classes
    # lty.up: The line type before the cut of the dendrogram. 
    # Line types can either be specified as an integer 
    # (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) 
    # lty.down: The line type for the classes (below the cut of the dendrogram)
    # (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) 
    # lwd.up: The line width beffore the cut of the dendrogram, 
    # Line width must be a positive number, defaulting to 1. 
    # lwd.down: The line width for the classes (below the cut of the dendrogram)
    # type: The type of dendrogram. Can be either "rectangle" or "triangle"
    # knot.pos: Possition of the knots; c("mean","bary","left","right","random")
    # show.label: logical value indicating whether the labels of the objects should be printed
    # only.tree: logical value indicating whether only the dendrogram should be printed
    # members: NULL or a vector with length size of a dissimilarity structure as produced by dist
    # ...   Further graphical arguments
 
    # ==================== Checking function arguments ===================
    if(missing(members)) members <- NULL
    opar <- par(no.readonly=TRUE)
    KNOTS <- c("mean","bary","left","right","random")
    knot.pos <- match.arg(knot.pos, KNOTS)
    TYPES <- c("rectangle", "triangle")
    type <- match.arg(type, TYPES)
    if (mode(k)!="numeric" || length(k)!=1 || k<=1 || (k%%1)!=0)
        stop("Invalid number of classes 'k'. Must be an integer larger than 1")
       
    ._a2r_counter    <<- 0
    ._a2r_hclu       <<- x
    ._a2r_envir      <<- environment()
    nn <- length(x$order) - 1
    ._a2r_height_cut <<- mean(x$height[nn-k+1:2])
    ._a2r_group      <<- 0    
    n.indiv <- length(x$order)
    groups.o <- .cutree.order(x, k=k)[x$order]    
    bottom <- if(is.null(members)) 0 else x$height[nn] * -.2     
    if (only.tree){
        if (is.null(members)) 
            plot(0, type="n", xlim=c(0.5,n.indiv+.5), ylim=c(bottom,x$height[nn]), 
                 xaxs="i", axes=FALSE, xlab="", ylab="") 
        else   plot(0, type="n", xlim=c(0.5,sum(members)+.5), ylim=c(bottom,x$height[nn]), 
                    xaxs="i", axes=FALSE, xlab="", ylab="")
        # call to the ** recursive function ** .rec.hclust
        .rec.hclust(nn, col=col.up, lty=lty.up, lwd=lwd.up)        
        axis(2)
        return(NULL)
    }    
    # prepare the layout
    matlayout <- matrix(c(1,2), nc=1, nr=2)
    widths    <- c(1,9)
    heights   <- c(9,1)
    layout(matlayout, width=widths, height=heights)        
    # Plotting the tree (1)
    par(mar=c(0,3,3,2))
    if(is.null(members)) plot(0,type="n",xlim=c(0.5,n.indiv+.5), ylim=c(bottom,x$height[nn]), xaxs="i", axes=FALSE, xlab="",ylab="") 
    else plot(0,type="n",xlim=c(0.5,sum(members)+.5), ylim=c(bottom,x$height[nn]), xaxs="i", axes=FALSE, xlab="",ylab="") 
    #call to the ** recursive function ** .rec.hclust
    .rec.hclust(nn, col=col.up, lty=lty.up, lwd=lwd.up)
    title(c("REBUS", paste("Colored Dendrogram of Outer and Inner Residuals (",k,"classes )")), cex.main=1)
    axis(2)
    # Plotting name of observations only if show.labels=TRUE (2)   
    if (show.labels){
        par(mar=c(2,2.5,0,2))
        par(srt=90)
        obs.labels <- toupper(substr(x$labels[x$order],1,6))
        if (is.null(members)) {
            plot(0,type="n", xlim=c(0.5,n.indiv+.5), ylim=c(0,1), xaxs="i", axes=FALSE, xlab="",ylab="") 
            text(1:n.indiv, 0, obs.labels, pos=4, col=col.down[groups.o], cex=.7)
        }
        else{
            plot(0,type="n",xlim=c(0.5,sum(members)+.5), ylim=c(0,1), xaxs="i", axes=FALSE, xlab="",ylab="") 
            xo <- members[x$order]
            text(cumsum(xo)-xo/2, 0, obs.labels, pos=4, col=col.down[groups.o])
        }
    }
    par(opar) # reset parameter
}

