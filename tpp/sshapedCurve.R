xc <- c(0.8,0.5,0.2,-0.2,-0.5,-0.8)
yn <- c(0.02,0.10,0.29,0.20,0.80,1.00)
nls.mod1 <- nls(yn~1/(1+exp(b*xc-a)),start=list(a=0,b=10))
summary(nls.mod1)

plot(xc,yn,xlim=c(-1,1))
x.grid <- seq(-1,1,by=0.01)
lines(x.grid,predict(nls.mod1,newdata=list(xc=x.grid)))
