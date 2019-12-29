rm(list = ls())
sigmoid = function(x) {
  1 / (1 + exp(-x))
}
sigmoid = function(x) {
   1 / (1+exp(-(1/(x))))
}
x <- seq(-5, 5, 0.01)
plot(x, sigmoid(x), col='blue')

x <- rnorm(100, 0, 5)
inv_logit <- function(x) {
  return(1 / (1 + exp(- x)))
}
y <- inv_logit(x)
plot(y ~ x)

sshape = function(x){
  1 / (1 + exp(x))
}
x <- seq(-5, 5, 0.01)
plot(x, sshape(x), col='red')


