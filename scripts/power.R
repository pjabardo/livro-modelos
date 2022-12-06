powerLaw <- function(p=2.5){

    x1 <- seq(0.1, 1, len=91)
    x2 <- seq(1, 10, len=91)
    x3 <- seq(10, 100, len=91)

    y1 <- x1^p
    y2 <- x2^p
    y3 <- x3^p

    pdf("../figuras/power.pdf", width=9, height=3.5)
    par(mfrow=c(1,3))
    plot(x1, y1, ty='l', xlab='x', ylab='y'); abline(h=max(y1))
    plot(x2, y2, ty='l', xlab='x', ylab='y'); abline(h=max(y2))
    plot(x3, y3, ty='l', xlab='x', ylab='y'); abline(h=max(y3))
    dev.off()
}
