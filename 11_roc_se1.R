roc_1 <- read.dta("data/roc_zm63_1.dta")
roc_2 <- read.dta("data/roc_zm63_2.dta")
roc <- rbind(roc_1, roc_2)


roc_p_1 <- roc(roc_1$ref, roc_1$prob)
roc_p_2 <- roc(roc_2$ref, roc_2$prob)
plot(roc_p_1, col=1, lty=c(2), main="", grid=c(0.2, 0.2))
plot(roc_p_2, col=2, lty=1, lwd=1, add=TRUE)
legend(x=0.1, y=0.85,unique(roc$survey),pch=15, col=1:length(unique(roc$survey)), bty="n", cex=0.6)
legend("bottomright",c("male","female"),lty=2:1, col=1:2, bty="n", cex=0.6)
