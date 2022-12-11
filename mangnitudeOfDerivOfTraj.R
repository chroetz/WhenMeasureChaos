u <- attractor$u[1e6:(1e6+4e4),]
plot(u[,1], type="l")
du <- apply(u, 2, diff)*1000
plot(du[,1], type="l")
ndu <- sqrt(rowMeans(du^2))
plot(ndu, type="l", ylim=c(0, 159))
grid()
avg1 <- double(length(ndu))
avg2 <- double(length(ndu))
avg3 <- double(length(ndu))
k1 <- length(ndu)/100
k2 <- length(ndu)/20
k3 <- length(ndu)/4
for (i in seq_along(avg1)) {
  avg1[i] <- mean(ndu[pmax(1, i-k1):pmin(length(ndu), i+k1)])
  avg2[i] <- mean(ndu[pmax(1, i-k2):pmin(length(ndu), i+k2)])
  avg3[i] <- mean(ndu[pmax(1, i-k3):pmin(length(ndu), i+k3)])
}
plot(ndu, type="l", ylim=c(0, max(ndu)))
grid()
lines(avg1, col=2, lwd=2)
lines(avg2, col=3, lwd=2)
lines(avg3, col=4, lwd=2)


plt <- ggplot(tibble(t = (1:4e4)/1e3, speed = ndu)) + geom_line(aes(x=t, y=speed)) + ylim(c(0, NA))
plt
ggsave("speed.pdf", plot=plt, width=8, height=4)
