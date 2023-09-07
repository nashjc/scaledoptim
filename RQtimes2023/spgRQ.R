# spgRQ.R
mbt<-2 # default replication in microbenchmark
require(microbenchmark)
molerfast <- function(n) {
  # A fast version of `molermat'
  A <- matrix(0, nrow = n, ncol = n)
  j <- 1:n
  for (i in 1:n) {
    A[i, 1:i] <- pmin(i, 1:i) - 2
  }
  A <- A + t(A)
  diag(A) <- 1:n
  A
}

axmolerfast <- function(x, AA=1) {
  # A fast and memory-saving version of A%*%x  
  # For Moler matrix. Note we need a matrix argument to match other functions
  n <- length(x)
  j <- 1:n
  ax <- rep(0, n)
  for (i in 1:n) {
    term <- x * (pmin(i, j) - 2)
    ax[i] <- sum(term[-i]) 
  }
  ax <- ax + j*x
  ax
}
rqfast<-function(x){
  rq<-as.numeric(t(x) %*% axmolerfast(x))
  rq
}
rqneg<-function(x) { -rqfast(x)}
proj <- function(x) {sign(x[1]) * x/sqrt(c(crossprod(x))) } # from ravi
# Note that the c() is needed in denominator to avoid error msgs
require(BB)
n<-100
x<-rep(1,n)
x<-x/as.numeric(sqrt(crossprod(x)))
AA<-molerfast(n)
teig<-microbenchmark(evs<-eigen(AA), times=mbt)$time
cat("eigen time =", mean(teig)*0.001,"\n")
tmin<-microbenchmark(amin<-spg(x, fn=rqfast, project=proj, 
                               control=list(trace=TRUE)), times=mbt)$time
#amin
tmax<-microbenchmark(amax<-spg(x, fn=rqneg, project=proj, 
                               control=list(trace=FALSE)), times=mbt)$time
#amax
evalmax<-evs$values[1]
evecmax<-evs$vectors[,1]
evecmax<-sign(evecmax[1])*evecmax/sqrt(as.numeric(crossprod(evecmax))) # normalize
emax<-list(evalmax=evalmax, evecmax=evecmax)
# save(emax, file="temax.Rdata")
evalmin<-evs$values[n]
evecmin<-evs$vectors[,n]
evecmin<-sign(evecmin[1])*evecmin/sqrt(as.numeric(crossprod(evecmin)))
avecmax<-amax$par
avecmin<-amin$par
avecmax<-sign(avecmax[1])*avecmax/sqrt(as.numeric(crossprod(avecmax)))
avecmin<-sign(avecmin[1])*avecmin/sqrt(as.numeric(crossprod(avecmin)))
cat("minimal eigensolution: Value=",amin$value,"in time ",mean(tmin)*0.001,"\n")
cat("Eigenvalue - result from eigen=",amin$value-evalmin,"  vector max(abs(diff))=",
    max(abs(avecmin-evecmin)),"\n\n")
#print(amin$par)
cat("maximal eigensolution: Value=",-amax$value,"in time ",mean(tmax)*0.001,"\n")
cat("Eigenvalue - result from eigen=",-amax$value-evalmax,"  vector max(abs(diff))=",
    max(abs(avecmax-evecmax)),"\n\n")
#print(amax$par)

# require(compiler)
nmax<-1
stable<-matrix(NA, nrow=nmax, ncol=4) # to hold results
# =========== works to here, but spg is slower than eigen
# loop over sizes
for (ni in 1:nmax){
  ni<-1
  n<-50*ni
  x<-runif(n) # generate a vector 
  AA<-molerfast(n) # make sure defined
  stable[[ni, 1]]<-n
  tbld<-microbenchmark(AA<-molerfast(n), times=mbt)
  tspg<-microbenchmark(aspg<-spg(x, fn=rqneg, project=proj, 
                                 control=list(trace=FALSE)), times=mbt)
  teig<-microbenchmark(aseig<-eigen(AA), times=mbt)
  stable[[ni, 2]]<-mean(tspg$time)*0.001
  stable[[ni, 3]]<-mean(tbld$time)*0.001
  stable[[ni, 4]]<-mean(teig$time)*0.001
}
spgtym<-data.frame(n=stable[,1], spgrqt=stable[,2], tbld=stable[,3], teig=stable[,4])
print(round(spgtym,0))
