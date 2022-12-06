koch01tmp <- function() list(x=c(0, 1/3, 1/2, 2/3, 1), y=c(0, 0, 1/(2*sqrt(3)), 0, 0))
koch01 <- function() list(x=c(1/3, 1/2, 2/3), y=c(0, 1/(2*sqrt(3)), 0))

applyKoch <- function(x, y, N=0, k=koch01()){
    
    dx <- x[2]-x[1]
    dy <- y[2]-y[1]

    l <- sqrt(dx*dx+dy*dy)
    ih <- c(dx/l, dy/l)
    jh <- c(-ih[2], ih[1])

    
    x2 <- c(x[1], x[1]+l*(ih[1]*k$x + jh[1]*k$y), x[2])
    y2 <- c(y[1], y[1]+l*(ih[2]*k$x + jh[2]*k$y), y[2])
    #lines(x2, y2, ty='l')
    if (N > 0){
        xx <- x2
        yy <- y2
        x0 <- double(0)
        y0 <- double(0)
        for (i in 1:4){
            pp <- applyKoch(c(xx[i], xx[i+1]), c(yy[i], yy[i+1]), N-1, k)
            np <- length(pp$x)
            #lines(pp$x, pp$y)
            x0 <- c(x0, pp$x[1:(np-1)])
            y0 <- c(y0, pp$y[1:(np-1)])
        }
        return(list(x=c(x0, x[2]), y=c(y0, y[2])))
    }
                       
    return(list(x=x2, y=y2))
}


kochFlake <- function(N=1, p=NULL){
    if (is.null(p)){
        x <- c(0, 0.5, 1, 0)
        y <- c(0, sqrt(3)/2, 0, 0)
        p <- list(x=x, y=y)
    }
    x <- p$x
    y <- p$y
    if (N==0) return(p)


    
    x2 <- double(0)
    y2 <- double(0)

    
    np <- length(p$x)
    k <- koch01()
    for (i in 1:(np-1)){
        pp <- applyKoch(c(x[i], x[i+1]), c(y[i], y[i+1]), N-1, k)
        n2 <- length(pp$x)-1
        x2 <- c(x2, pp$x[1:n2])
        y2 <- c(y2, pp$y[1:n2])
    }
    
    return(list(x=c(x2, x[1]), y=c(y2, y[1])))
}
    
mylen <- function(p){
    x <- p$x
    y <- p$y

    np <- length(x)

    dx <- x[2:np] - x[1:(np-1)]
    dy <- y[2:np] - y[1:(np-1)]

    return(sum(sqrt(dx*dx+dy*dy)))
    
}
myarea <- function(p){
    
    x <- p$x
    y <- p$y

    np <- length(x)
    return(0.5*abs(sum(x[1:(np-1)]*y[2:np] - y[1:(np-1)]*x[2:np])))
}

kochFigure <- function(){

    p <- lapply(0:5, kochFlake)

    xlim <- range(p[[6]]$x)
    ylim <- range(p[[6]]$y)
    pdf("../figuras/koch.pdf", width=9, height=3.5)
    par(mfrow=c(2,3), mar=rep(0.5, 4))
    for (i in 1:6){
        plot.new()
        plot.window(xlim, ylim, asp=1)
        lines(p[[i]]$x, p[[i]]$y)
        legend("topleft", legend=paste("N =", i-1), bty='n', cex=2)
    }
    dev.off()
}


selfSimilarity <- function(p, N=9){
    if (is.null(p)) p <- kochFlake(N)

    xlims <- matrix(c(-0.0737166, -0.04600082, -0.01116245, -0.002867262, -0.0004303842,
                      0.2865885,  0.10013692,  0.03401645,  0.009862324, 0.0012384193), nr=2, nc=5, byrow=TRUE)
    ylims <- matrix(c( 0.3187060, 0.4801025, 0.5479675, 0.5691191, 0.5761624, 
                      0.6919354, 0.6288899, 0.5895606, 0.5822905, 0.5780322), nr=2, nc=5, byrow=TRUE)
    x <- p$x
    y <- p$y
    xlims=cbind(range(x), xlims)
    ylims=cbind(range(y), ylims)
    png('../figuras/koch-self.png', width=1920, height=1280, pointsize=30)
    par(mfrow=c(2,3), mar=c(4,4,1,1))
    for (i in 1:6){
        
        plot(p$x, p$y, ty='l', asp=1, xlab='x', ylab='y', xlim=xlims[,i], ylim=ylims[,i], lwd=2)
        if (i < 6) rect(xlims[1,i+1], ylims[1, i+1], xlims[2,i+1], ylims[2,i+1], lwd=2)
    }
    dev.off()
}
