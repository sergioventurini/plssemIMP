n <- 1000
k <- 1:999
p <- k/n
m <- n*p
sd <- sqrt(n*p*(1 - p))
L <- 2
a <- round(m - L*sd, 0)
b <- round(m + L*sd, 0)
# a <- m - L*sd
# b <- m + L*sd
# probs <- 1 - (pbinom(b, n, p) - pbinom(a, n, p))
probs <- 1 - pbinom(b, n, p)
B <- 1000
plot(k, probs*B, type = "l")
k <- 500
p <- k/n
m <- n*p
sd <- sqrt(n*p*(1 - p))
a <- round(m - L*sd, 0)
b <- round(m + L*sd, 0)
(1 - pbinom(b - 1, n, p))*B
abline(v = k, lty = 2)
abline(h = (1 - pbinom(b, n, p))*B, lty = 2)
