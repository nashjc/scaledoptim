knitr::opts_chunk$set(echo = TRUE)
## require(bookdown) # language engine to display text - does not seem necessary

require(optimx, quietly=TRUE)
pr <- function(y) {
- prod(y)*(1-sum(y))
}
cat("test the simple product for n=5\n")
meth <- c("Nelder-Mead", "BFGS")
n<-5
  st<-1:(n-1)/(n*n)
   ans<-opm(st, pr, gr="grcentral", control=list(trace=0))
   ao<-summary(ans,order=value)
print(ao)

nll <- function(y) {
  if ((any(y <= 10*.Machine$double.xmin)) || (sum(y)>1-.Machine$double.eps))
         .Machine$double.xmax
  else   - sum(log(y)) - log(1-sum(y))
}
nll.g <- function(y) { - 1/y + 1/(1-sum(y))} # so far not safeguarded
n<-5
x0<-(2:n)/n^2
library(numDeriv)
dx0<-nll.g(x0)
dx0n<-grad(nll, x0)
cat("Max Abs diff between analytic and approx. gradient =",max(abs(dx0-dx0n)),"\n")

require(optimx, quietly=TRUE)
mset<-c("L-BFGS-B", "BFGS", "CG", "spg", "ucminf", "nlm", "nlminb", "nvm", "ncg", "tnewt")
# numerical approximation using forward difference
a5<-opm(2:n/n^2, nll, gr="grfwd", method=mset, 
          control=list(dowarn=FALSE, kkt=FALSE,trace=0))
summary(a5, order=value)
# analytical gradient
a5g<-opm(2:n/n^2, nll, nll.g, method=mset, 
          control=list(dowarn=FALSE,kkt=FALSE, trace=0))
summary(a5g, order=value)
# analytical gradient and bounds on parameters 
a5gb<-opm(2:n/n^2, nll, nll.g, lower=0, upper=1, method=mset, 
          control=list(dowarn=FALSE,kkt=FALSE,trace=0))
summary(a5gb, order=value)

require(BB, quietly=TRUE)
nllrv <- function(x) {- sum(log(x))}
nllrv.g <- function(x) {- 1/x }
proj <- function(x) {x/sum(x)}
n <- 5
aspg <- spg(par=(1:n)/n^2, fn=nllrv, gr=nllrv.g, project=proj, 
            control=list(trace=TRUE, triter=1))
aspgn <- spg(par=(1:n)/n^2, fn=nllrv, project=proj, 
             control=list(trace=TRUE, triter=1)) # using internal grad approx.
cat("F_optimal: with gradient=",aspg$value,"  num. approx.=",aspgn$value,"\n")
pbest<-rep(1/n, n)
cat("fbest = ",nllrv(pbest),"  when all parameters = ", pbest[1],"\n")
cat("deviations:  with gradient=",max(abs(aspg$par-pbest)),
    "   num. approx.=",max(abs(aspg$par-pbest)),"\n")

enll <- function(lx) {
    x<-exp(lx)
    fval<-  - sum( log( x/sum(x) ) ) 
}
enll.g <- function(lx){
    x<-exp(lx)
    g<-length(x)/sum(x) - 1/x
    gval<-g*exp(lx)
}

require(optimx, quietly=TRUE) # just to be sure
st<-1:5/10 # 5 parameters, crude scaling to start
a5x<-opm(st, enll, enll.g, method="MOST", control=list(trace=0))
a5x<-a5x[order(a5x$value),]
cat("Proposed best solution has minimum=",a5x[1,length(st)+1],"\n")
cat("Coeffs:")
print(a5x[1,1:length(st)])

require(optimx, quietly=TRUE)
st<-1:100/1e3 # large
stenll<-enll(st)
cat("Initial function value =",stenll,"\n")
tym<-system.time(acgbig<-Rcgmin(st, enll, enll.g, 
                                control=list(trace=0, tol=1e-32)))[[3]]
cat("Time = ",tym,"  fval=",acgbig$value,"\n")
xnor<-acgbig$par/sum(acgbig$par)
print(xnor)

library(optimx)
proj2 <- function(theta) {
    theta2 <- theta^2
    s2 <- theta2 / (1 + theta2)
    cumprod(c(1, s2)) * c(1-s2, 1)
 }
obj <- function(theta) - sum(log(proj2(theta)))
 n <- 5
 ans <- spg(seq(n-1), obj)
 proj2(ans$par)

n<-100
# check 
obj(seq(n-1))
ans100 <- spg(seq(n-1), obj, control=list(trace=FALSE), quiet=TRUE)
ans100$value
proj2( (ans100$par) )

# turn off kkt to save time
tmeth<-c("ncg", "nvm", "lbfgs", "ucminf", "bobyqa", "tnewt", "slsqp")
mfv<-10*(n-1)^2 # set fn eval count big enough to avoid commonArgs error (bobyqa)
allansf<- opm(seq(n-1)/n, obj, gr="grfwd", method=tmeth, 
              control=list(kkt=FALSE, dowarn=FALSE,maxfeval=mfv))
summary(allansf, order = "list(round(value, 3), fevals)", par.select = FALSE)
allansc<- opm(seq(n-1)/n, obj, gr="grcentral", method=tmeth, 
              control=list(kkt=FALSE, dowarn=FALSE,maxfeval=mfv))
summary(allansc, order = "list(round(value, 3), fevals)", par.select = FALSE)
allansp<- opm(seq(n-1)/n, obj, gr="grpracma", method=tmeth, 
              control=list(kkt=FALSE, dowarn=FALSE,maxfeval=mfv))
summary(allansp, order = "list(round(value, 3), fevals)", par.select = FALSE)

molerbuild<-function(n){ # Create the moler matrix of order n
   # A[i,j] = i for i=j, min(i,j)-2 otherwise
   A <- matrix(0, nrow = n, ncol = n)
   j <- 1:n
   for (i in 1:n) {
      A[i, 1:i] <- pmin(i, 1:i) - 2
   }
   A <- A + t(A)
   diag(A) <- 1:n
   A
}

raynum<-function(x, A){
   rayquo<-as.numeric((t(x)%*%A)%*%x)
}

proj<-function(x) { x/sqrt(sum(x*x)) }

require(BB, quietly=TRUE)
n<-10
set.seed(4321)
x<-runif(n) 
B<-molerbuild(n)
tmin<-system.time(asprqmin<-spg(x, fn=raynum, project=proj, A=B,
    control=list(trace=TRUE, triter=1,maxit=1000)))[[3]]
tmax<-system.time(asprqmax<-spg(x, fn=raynum, project=proj, A=-B,
    control=list(trace=TRUE, triter=1,maxit=1000)))[[3]]
xy <- rep(1/sqrt(n), n)
tmax2<-system.time(asprqmax2<-spg(xy, fn=raynum, project=proj, A=-B,
    control=list(trace=TRUE, triter=1,maxit=1000)))[[3]]
cat("maximal eigensolution: Value=",-asprqmax2$value,"in time ",tmax2,"\n")
print(asprqmax2$par)
cat("minimal eigensolution: Value=",asprqmin$value,"in time ",tmin,"\n")
print(asprqmin$par)
eigen(B)$values

ssums<-function(x){
  n<-length(x)
  tt<-sum(x)
  ss<-1:n
  xx<-(x/tt)*(x/tt)
  sum(ss*xx)
}

cat("Try penalized sum\n")
require(optimx)
st<-runif(3)
aos<-opm(st, ssums, gr="grcentral", method="MOST")
# rescale the parameters
nsol<-dim(aos)[1]
for (i in 1:nsol){ 
  tpar<-aos[i,1:3] 
  ntpar<-sum(tpar)
  tpar<-tpar/ntpar
#  cat("Method ",aos[i, "meth"]," gives fval =", ssums(tpar))
  aos[i, 1:3]<-tpar 
}
summary(aos,order=value)[1:5,]

ssum<-function(x){
  n<-length(x)
  ss<-1:n
  xx<-x*x
  sum(ss*xx)
}
proj.simplex <- function(y) {
# project an n-dim vector y to the simplex Dn
# Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}
# Ravi Varadhan, Johns Hopkins University
# August 8, 2012
n <- length(y)
sy <- sort(y, decreasing=TRUE)
csy <- cumsum(sy)
rho <- max(which(sy > (csy - 1)/(1:n)))
theta <- (csy[rho] - 1) / rho
return(pmax(0, y - theta))
}
as<-spg(st, ssum, project=proj.simplex)
cat("Using project.simplex with spg: fmin=",as$value," at \n")
print(as$par)

library(quadprog)
Dmat<-diag(c(1,2,3))
Amat<-matrix(c(1, 1, 1), ncol=1)
bvec<-c(1)
meq=1
dvec<-c(0, 0, 0)
ans<-solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)
ans
