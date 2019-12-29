load(file = file.path(resDir,'output_step1.RData'))
load(file = file.path(resDir,'output_step2.RData'))
temps <- as.numeric(substring(tempRange,5))
ratios <- unname(unlist(df[1,5:14]))
x <- temps
y <- ratios
fit <- list()
st.coarse <- expand.grid(p=c(0,0.3),k=seq(0,4000,by=1000),m=seq(30,60,by=15))
st.fine   <- expand.grid(p=c(0,0.3),k=seq(0,8000,by=200),m=seq(30,80,by=10))
for (st in list(st.coarse,st.fine)) {
  tryCatch( {
    mod <- nls2(y~sigmoid(p,k,m,x),data=list(x=x,y=y),start=st,algorithm="brute-force",control=nls.control(warnOnly=T,maxiter=5000))
    fit <-
      nls2(y~sigmoid(p,k,m,x),data=list(x=x,y=y),start=mod,control=nls.control(warnOnly=F),algorithm="port",lower=c(0,1,10),upper=c(0.4,100000,100))
    obj <- list()
    obj$plat  <- as.numeric(coefficients(fit)[1])
    obj$k     <- as.numeric(coefficients(fit)[2])
    obj$tm    <- as.numeric(coefficients(fit)[3])
    obj$slope <- as.numeric(sigmoid.d1(obj$plat,obj$k,obj$tm,obj$tm))
    y.fit <- sigmoid(obj$plat,obj$k,obj$tm,temps)
    obj$y.fit <- y.fit
    obj$resid <- ratios - y.fit
    obj$r2 <- 1-(sum(obj$resid^2)/(length(ratios)*var(ratios)))
    obj$rmsd <- sqrt( sum(obj$resid^2)/length(ratios) )
    return(obj)
  },error = function(e) {})
}
