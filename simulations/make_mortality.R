f <- function (t) {
    m <- 0.5 + (0.95-0.5) * (dnorm(tt, mean=160, sd=10) + dnorm(tt, mean=260, sd=10) + dnorm(tt, mean=210, sd=30)*0.8) / (dnorm(0, sd=10))
    return(m)
}

tt <- 1:365

pdf(file="survival.pdf", width=6, height=4, pointsize=10)
plot(tt, f(tt), type='l', xlab='day of the year', ylab='probability of survival', ylim=c(0,1))
abline(h=c(0.5, 0.95), lty=3)
dev.off()


out <- data.frame(t=tt, p=f(tt))
write.csv(out, file="survival.csv", row.names=FALSE)

