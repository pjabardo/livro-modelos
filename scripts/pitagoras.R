pitagoras <- function(){


    a = 5
    phi <- pi/6
    b <- a*cos(phi)
    h <- b*sin(phi)
    x <- b*cos(phi)
    cc <- sqrt(a*a-b*b)
    pdf("../figuras/pitagoras.pdf", width=7, height=3)
    plot.new()
    par(mar=c(0,0,0,0))
    plot.window(c(0,a), c(-0.5,h*1.1), asp=1)

    polygon(c(0, a, x), c(0, 0, h), lwd=2)

    lines(c(x,x), c(h,0), lwd=1)
    text(a/2, 0, "a", pos=1, cex=2)
    text(b/2*cos(phi), b/2*sin(phi)+0.05, "b", pos=3, cex=2)
    text(a-cc/2*cos(pi/2-phi), cc/2*sin(pi/2-phi), "c", pos=4, cex=2)
    th <- seq(0, phi, len=11)
    r <- 0.5
    lines(r*cos(th), r*sin(th))
    text(r*cos(phi/2), r*sin(phi/2), expression(phi), pos=4, cex=1.7)
    lines(r*sin(th)+x, -r*cos(th)+h)
    text(r*sin(phi/2)+x, -r*cos(phi/2)+h, expression(phi), pos=1, cex=1.5)
    text(mean(c(x,x,0)), mean(c(0,0,h)), expression(A[b]), cex=2)
    text(mean(c(x,x,a)), mean(c(0,0,h)), expression(A[c]), cex=1.5)
    text(a/2, h, expression(A[a]), cex=2.5)
    dev.off()
}

