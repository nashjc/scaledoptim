## ----setup, echo=FALSE---------------------------------------------------------------------------------------------------------------------------
mbt<-25 # default to 5 repetitions in microbenchmark while sorting out text
msect<-function(times){
   round(mean(times)*0.001,0)
}
msecr<-function(times){
#   round((max(times)-min(times))*0.001,0)
   round(sd(times)*0.001,0)
}


## ----molermat, echo=TRUE-------------------------------------------------------------------------------------------------------------------------
molermat<-function(n){
   A<-matrix(NA, nrow=n, ncol=n)
   for (i in 1:n){
      for (j in 1:n) {
          if (i == j) A[i,i]<-i
          else A[i,j]<-min(i,j) - 2
      }
   }
   A
}


## ----molerfast, echo=TRUE------------------------------------------------------------------------------------------------------------------------
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


## ----mnolertime, echo=FALSE, cache=TRUE----------------------------------------------------------------------------------------------------------
require(microbenchmark)
nmax<-5
mtable<-matrix(NA, nrow=nmax, ncol=5) # to hold results
rtable<-mtable
# loop over sizes
for (ni in 1:nmax){
  n<-100*ni
  mtable[[ni, 1]]<-n
  rtable[[ni, 1]]<-n
  # Note "unit" argument is ONLY for display. time is in nanoseconds
  ti<-microbenchmark(ai<-molermat(n), unit="us", times=mbt)$time
  tfi<-microbenchmark(afi<-molerfast(n), unit="us", times=mbt)$time
  if (! identical(ai, afi)) stop("Different outcomes == molermat, molerfast")
  matsize<-as.numeric(object.size(ai))
  tevs<-microbenchmark(evs<-eigen(ai), unit="us", times=mbt)$time
  mtable[[ni,2]]<-msect(ti) 
  mtable[[ni,3]]<-matsize
  mtable[[ni,4]]<-msect(tevs)
  mtable[[ni,5]]<-msect(tfi)
  rtable[[ni,2]]<-msecr(ti) 
  rtable[[ni,3]]<-matsize
  rtable[[ni,4]]<-msecr(tevs)
  rtable[[ni,5]]<-msecr(tfi)
}


## ----molertimedisp, echo=FALSE-------------------------------------------------------------------------------------------------------------------
bmattym<-data.frame(n=mtable[,1], matsize=mtable[,3], buildi=mtable[,2],
     buildir=rtable[,2], eigentime=mtable[,4], eigentimr=rtable[,4], 
     bfast=mtable[,5], bfastr=rtable[,5])
print(bmattym)
cat("matsize - matrix size in bytes\n")
cat("eigentime - all eigensolutions time\n")
cat("buildi - interpreted build time, range\n")
cat("bfast - interpreted vectorized build time\n")
cat("Times converted to milliseconds\n")


## ----drawtime1, echo=FALSE-----------------------------------------------------------------------------------------------------------------------
ti<-as.vector(mtable[,2])
tf<-as.vector(mtable[,5])
matsize<-as.vector(mtable[,3])
n<-as.vector(mtable[,1])
plot(n, ti)
xx<-1:max(mtable[,1])
n2<-n*n
itime<-lm(ti~n+n2)
summary(itime)
ftime<-lm(tf~n+n2)
summary(ftime)
iti<-coef(itime)
yy<-iti[1]+iti[2]*xx+iti[3]*xx*xx
points(xx,yy, type='l')
fti<-coef(ftime)
ww<-fti[1]+fti[2]*xx+fti[3]*xx*xx
points(n, tf, col='red')
points(xx, ww, type='l', col='red')
title(main="Execution time vs matrix size")
title(sub="molermat (black) and molerfast (red) matrix builds")


## ----drawtime2, echo=FALSE-----------------------------------------------------------------------------------------------------------------------
matsizmod<-lm(matsize~n+n2)
summary(matsizmod)
cos<-coef(matsizmod)
zz<-cos[1]+cos[2]*xx+cos[3]*xx*xx
plot(n, matsize)
points(xx, zz, type='l')
title(main="Execution time vs matrix size")
title(sub="eigen() on Moler matrix")


## ----rqdir, echo=TRUE----------------------------------------------------------------------------------------------------------------------------
rqdir<-function(x, AA){
  rq<-0.0
  n<-length(x) # assume x, AA conformable
  for (i in 1:n) {
     for (j in 1:n) {
        rq<-rq+x[i]*AA[[i,j]]*x[j]
     }
  }
  rq
}


## ----raynum1, echo=TRUE--------------------------------------------------------------------------------------------------------------------------
ray1<-function(x, AA){
    rq<-  t(x)%*%AA%*%x
}


## ----raynum2, echo=TRUE--------------------------------------------------------------------------------------------------------------------------
ray2<-function(x, AA){
    rq<-  as.numeric(crossprod(x, crossprod(AA,x)))
}


## ----raynum3, echo=TRUE--------------------------------------------------------------------------------------------------------------------------
ray3<-function(x, AA, ax=axftn){
    # ax is a function to form AA%*%x 
    rq<- - as.numeric(crossprod(x, ax(x, AA)))
}


## ----axm, echo=TRUE------------------------------------------------------------------------------------------------------------------------------
ax<-function(x, AA){
   u<- as.numeric(AA%*%x)
}

axx<-function(x, AA){
   u<- as.numeric(crossprod(AA, x))
}


## ----aximp, echo=TRUE----------------------------------------------------------------------------------------------------------------------------
aximp<-function(x, AA=1){ # implicit moler A*x
   n<-length(x)
   y<-rep(0,n)
   for (i in 1:n){
      tt<-0.
      for (j in 1:n) {
          if (i == j) tt<-tt+i*x[i]
          else tt<-tt+(min(i,j) - 2)*x[j]
      }
      y[i]<-tt 
   }
   y
}
ident<-function(x, B=1) x # identity


## ----axmfcode, echo=TRUE-------------------------------------------------------------------------------------------------------------------------
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


## ----axftn, echo=TRUE----------------------------------------------------------------------------------------------------------------------------
dyn.load("moler.so")
cat("Is the mat multiply loaded? ",is.loaded("moler"),"\n")

axftn<-function(x, AA=1) { # ignore second argument
   n<-length(x) # could speed up by having this passed
   vout<-rep(0,n) # purely for storage
   res<-(.Fortran("moler", n=as.integer(n), x=as.double(x), vout=as.double(vout)))$vout
}


## ----timeax1, echo=TRUE, cache=TRUE--------------------------------------------------------------------------------------------------------------
dyn.load("moler.so")
cat("Is the mat multiply loaded? ",is.loaded("moler"),"\n")
require(microbenchmark)
nmax<-5
ptable<-matrix(NA, nrow=nmax, ncol=11) # to hold results
# loop over sizes
for (ni in 1:nmax){
  n<-100*ni
  x<-runif(n) # generate a vector 
  ptable[[ni, 1]]<-n
  AA<-molermat(n)
  tax<- microbenchmark(oax<-ax(x, AA), times=mbt)$time
  taxx<-microbenchmark(oaxx<-axx(x, AA), times=mbt)$time
  if (! identical(oax, oaxx)) stop("oaxx NOT correct")
  taxftn<-microbenchmark(oaxftn<-axftn(x, AA=1), times=mbt)$time
  if (! identical(oax, oaxftn)) stop("oaxftn NOT correct")
  taximp<-microbenchmark(oaximp<-aximp(x, AA=1), times=mbt)$time
  if (! identical(oax, oaximp)) stop("oaximp NOT correct")
  taxmfi<-microbenchmark(oaxmfi<-axmolerfast(x, AA=1), times=mbt)$time
  if (! identical(oax, oaxmfi)) stop("oaxmfi NOT correct")
  ptable[[ni, 2]]<-msect(tax); ptable[[ni,3]]<-msecr(tax)
  ptable[[ni, 4]]<-msect(taxx); ptable[[ni, 5]]<-msecr(taxx)
  ptable[[ni, 6]]<-msect(taxftn); ptable[[ni, 7]]<-msecr(taxftn)
  ptable[[ni, 8]]<-msect(taximp); ptable[[ni,9]]<-msecr(taximp)
  ptable[[ni, 10]]<-msect(taxmfi); ptable[[ni,11]]<-msecr(taxmfi)
}

axtym<-data.frame(n=ptable[,1], ax=ptable[,2], sd_ax=ptable[,3],  axx=ptable[,4],
                  sd_axx=ptable[,5],  axftn=ptable[,6], sd_axftn=ptable[,7], 
                  aximp=ptable[,8], sd_aximp=ptable[,9], 
                  axmfast=ptable[,10], sd_axmfast=ptable[,11])
print(axtym)


## ----extabl1, echo=FALSE-------------------------------------------------------------------------------------------------------------------------
# explain table
expln <- c("ax = R matrix * vector  A %*% x",
   "axx = R crossprod A, x",
   "axftn = Fortran version of implicit Moler A * x",
   "aximp = implicit moler A*x in R",
   "axmfast = A fast and memory-saving version of A %*% x",
   "Times in milliseconds from microbenchmark")
for (exx in expln) { cat(exx,"\n")}


## ----adjaxtime, echo=FALSE-----------------------------------------------------------------------------------------------------------------------
cat("Times (in millisecs) adjusted for matrix build\n")
adjtym<-data.frame(n=axtym$n, axbld=axtym$ax+bmattym$buildi, 
     axxbld=axtym$axx+bmattym$buildi, 
     axftn=axtym$axftn, aximp=axtym$aximp)
print(adjtym)


## ----rqtime1, echo=FALSE, cache=TRUE-------------------------------------------------------------------------------------------------------------
dyn.load("moler.so")
  n<-500
  x<-runif(n) # generate a vector 
  AA<-molermat(n)
  tdi<-microbenchmark(rdi<-rqdir(x, AA))$time
  cat("Direct algorithm: ",msect(tdi),"sd=",msecr(tdi),"\n")
  t1i<-microbenchmark(r1i<-ray1(x, AA))$time
  cat("ray1: mat-mult algorithm: ",msect(t1i),"sd=",msecr(t1i),"\n")
  t2i<-microbenchmark(r2i<-ray2(x, AA))$time
  cat("ray2: crossprod algorithm: ",msect(t2i),"sd=",msecr(t2i),"\n")
  t3fi<-microbenchmark(r3i<-ray3(x, AA, ax=axftn))$time
  cat("ray3: ax Fortran + crossprod: ",mean(t3fi)*0.001,"\n")
  t3ri<-microbenchmark(r3i<-ray3(x, AA, ax=axmolerfast))$time
  cat("ray3: ax fast R implicit + crossprod: ",msect(t3ri),"sd=",msecr(t3ri),"\n")


## ----rayspg1, echo=TRUE, cache=TRUE--------------------------------------------------------------------------------------------------------------
# spgRQ.R
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
cat("eigen time =", msect(teig),"sd=",msecr(teig),"\n")
tmin<-microbenchmark(amin<-spg(x, fn=rqfast, project=proj, 
                               control=list(trace=FALSE)), times=mbt)$time
tmax<-microbenchmark(amax<-spg(x, fn=rqneg, project=proj, 
                               control=list(trace=FALSE)), times=mbt)$time
evalmax<-evs$values[1]
evecmax<-evs$vectors[,1]
evecmax<-sign(evecmax[1])*evecmax/sqrt(as.numeric(crossprod(evecmax))) # normalize
emax<-list(evalmax=evalmax, evecmax=evecmax)
evalmin<-evs$values[n]
evecmin<-evs$vectors[,n]
evecmin<-sign(evecmin[1])*evecmin/sqrt(as.numeric(crossprod(evecmin)))
emin<-list(evalmin=evalmin, evecmin=evecmin)
avecmax<-amax$par
avecmin<-amin$par
avecmax<-sign(avecmax[1])*avecmax/sqrt(as.numeric(crossprod(avecmax)))
avecmin<-sign(avecmin[1])*avecmin/sqrt(as.numeric(crossprod(avecmin)))
cat("minimal eigensolution: Value=",amin$value,"in time ",
      msect(tmin),"sd=",msecr(tmin),"\n")
cat("Eigenvalue - result from eigen=",amin$value-evalmin,"  vector max(abs(diff))=",
    max(abs(avecmin-evecmin)),"\n")
#print(amin$par)
cat("maximal eigensolution: Value=",-amax$value,"in time ",
     msect(tmax),"sd=",msecr(tmax),"\n")
cat("Eigenvalue - result from eigen=",-amax$value-evalmax,"  vector max(abs(diff))=",
    max(abs(avecmax-evecmax)),"\n")

# nmax<-5
# stable<-matrix(NA, nrow=nmax, ncol=4) # to hold results
# # =========== works to here, but spg is slower than eigen
# # loop over sizes
# for (ni in 1:nmax){
#   n<-50*ni
#   x<-runif(n) # generate a vector 
#   AA<-molerfast(n) # make sure defined
#   stable[[ni, 1]]<-n
#   tbld<-microbenchmark(AA<-molerfast(n), times=mbt)
#   tspg<-microbenchmark(aspg<-spg(x, fn=rqneg, project=proj, 
#                                  control=list(trace=FALSE)), times=mbt)
#   teig<-microbenchmark(aseig<-eigen(AA), times=mbt)
#   stable[[ni, 2]]<-msect(tspg$time)
#   stable[[ni, 3]]<-msect(tbld$time)
#   stable[[ni, 4]]<-msect(teig$time)
# }
# spgtym<-data.frame(n=stable[,1], spgrqt=stable[,2], tbld=stable[,3], teig=stable[,4])
# print(round(spgtym,0))


## ----runopx1, echo=TRUE, cache=TRUE--------------------------------------------------------------------------------------------------------------
require(optimx)
nobj<-function(x, AA=-AA){
   y<-x/sqrt(as.numeric(crossprod(x)))
   rq<- as.numeric(crossprod(y, crossprod(AA,y)))
}

ngrobj<-function(x, AA=-AA){
   y<-x/sqrt(as.numeric(crossprod(x))) 
   n<-length(x)
   dd<-sqrt(as.numeric(crossprod(x)))
   T1<-diag(rep(1,n))/dd
   T2<- x%o%x/(dd*dd*dd)
   gt<-T1-T2
   gy<- as.vector(2.*crossprod(AA,y))
   gg<-as.numeric(crossprod(gy, gt))
} 
mset<-c("L-BFGS-B", "BFGS", "ncg", "spg", "ucminf", "nlm", "nlminb", "nvm")
for (ni in 1:nmax){
  n<-20*ni
  x<-runif(n) # generate a vector 
  AA<-molerfast(n) # make sure defined
  aall<-opm(x, fn=nobj, gr=ngrobj, method=mset, AA=-AA, 
     control=list(trace=0,starttests=FALSE, dowarn=FALSE, kkt=FALSE))
  # optansout(aall, NULL)
  summary(aall, order=value, )
}


## ----rcgrun1, echo=TRUE,cache=TRUE---------------------------------------------------------------------------------------------------------------
ctable<-matrix(NA, nrow=10, ncol=2)
nmax<-5
for (ni in 1:nmax){
  n<-50*ni
  x<-runif(n) # generate a vector 
  AA<-molerfast(n) # define matrix
  tcgu<-microbenchmark(arcgu<-optimr(x, fn=nobj, gr=ngrobj, method="ncg",
          AA=-AA), times=mbt)
  ctable[[ni,1]]<-n
  ctable[[ni,2]]<-mean(tcgu$time)*0.001
}
cgtime<-data.frame(n=ctable[,1], tcgmin=ctable[,2])
print(round(cgtime,0))


## ----geradincode, echo=FALSE,cache=TRUE----------------------------------------------------------------------------------------------------------
ax<-function(x, AA){
   u<-as.numeric(AA%*%x)
}

bx<-function(x, BB){
   v<-as.numeric(BB%*%x)
}

geradin<-function(x, ax, bx, AA, BB, control=list(trace=TRUE, maxit=1000)){
# Geradin minimize Rayleigh Quotient, Nash CMN Alg 25
#  print(control)
  trace<-control$trace
  n<-length(x)
  tol<-n*n*.Machine$double.eps^2
  offset<-1e+5 # equality check offset
  if (trace) cat("geradin.R, using tol=",tol,"\n")
  ipr<-0 # counter for matrix mults
  pa<-.Machine$double.xmax
  R<-pa
  msg<-"no msg"
# step 1 -- main loop
  keepgoing<-TRUE
  while (keepgoing) {
    avec<-ax(x, AA); bvec<-bx(x, BB); ipr<-ipr+1
    xax<-as.numeric(crossprod(x, avec));  
    xbx<-as.numeric(crossprod(x, bvec));
    if (xbx <= tol) {
       keepgoing<-FALSE # not really needed
       msg<-"avoid division by 0 as xbx too small"
       break
    } 
    p0<-xax/xbx
    if (p0>pa) {
       keepgoing<-FALSE # not really needed
       msg<-"Rayleigh Quotient increased in step"
       break
    } 
    pa<-p0
    g<-2*(avec-p0*bvec)/xbx
    gg<-as.numeric(crossprod(g)) # step 6
    if (trace) cat("Before loop: RQ=",p0," after ",ipr," products, gg=",gg,"\n")
    if (gg<tol) { # step 7
       keepgoing<-FALSE # not really needed
       msg<-"Small gradient -- done"
       break
    } 
    t<- -g # step 8
    for (itn in 1:n) { # major loop step 9
       y<-ax(t, AA); z<-bx(t, BB); ipr<-ipr+1 # step 10
       tat<-as.numeric(crossprod(t, y)) # step 11
       xat<-as.numeric(crossprod(x, y)) 
       xbt<-as.numeric(crossprod(x, z)) 
       tbt<-as.numeric(crossprod(t, z)) 
       u<-tat*xbt-xat*tbt
       v<-tat*xbx-xax*tbt
       w<-xat*xbx-xax*xbt
       d<-v*v-4*u*w
       if (d<0) stop("Geradin: imaginary roots not possible") # step 13
       d<-sqrt(d) # step 14
       if (v>0) k<--2*w/(v+d) else k<-0.5*(d-v)/u
       xlast<-x # NOT as in CNM -- can be avoided with loop
       avec<-avec+k*y; bvec<-bvec+k*z # step 15, update
       x<-x+k*t
       xax<-xax+as.numeric(crossprod(x,avec))      
       xbx<-xbx+as.numeric(crossprod(x,bvec))      
       if (xbx<tol) stop("Geradin: xbx has become too small")
       chcount<-n - length(which((xlast+offset)==(x+offset)))
       if (trace) cat("Number of changed components = ",chcount,"\n")
       pn<-xax/xbx # step 17 different order
       if (chcount==0) {
         keepgoing<-FALSE # not really needed
         msg<-"Unchanged parameters -- done"
         break
       }
       if (pn >= p0) {
         if (trace) cat("RQ not reduced, restart\n")
         break # out of itn loop, not while loop (TEST!)
       }
       p0<-pn # step 19
       g<-2*(avec-pn*bvec)/xbx
       gg<-as.numeric(crossprod(g))
       if (trace) cat("Itn", itn," RQ=",p0," after ",ipr," products, gg=",gg,"\n")
       if (gg<tol){ # step 20
         if (trace) cat("Small gradient in iteration, restart\n")
         break # out of itn loop, not while loop (TEST!)
       }
       xbt<-as.numeric(crossprod(x,z)) # step 21
       w<-y-pn*z # step 22
       tabt<-as.numeric(crossprod(t,w))
       beta<-as.numeric(crossprod(g,(w-xbt*g)))
       beta<-beta/tabt # step 23
       t<-beta*t-g
    } # end loop on itn -- step 24
  } # end main loop -- step 25
# step 26
  ans<-list(x=x, RQ=p0, ipr=ipr, msg=msg)
}


## ----rungeradin10, echo=TRUE, cache=TRUE---------------------------------------------------------------------------------------------------------
cat("Test geradin with explicit matrix multiplication\n")
n<-10
AA<-molerfast(n)
BB=diag(rep(1,n))
x<-runif(n)
tg<-microbenchmark(ag<-geradin(x, ax, bx, AA=AA, BB=BB, 
   control=list(trace=FALSE)), times=mbt)
cat("Minimal eigensolution\n")
print(ag)
cat("Geradin time=",msect(tg$time),"sd=",msecr(tg$time),"\n")
tgn<-microbenchmark(agn<-geradin(x, ax, bx, AA=-AA, BB=BB,
   control=list(trace=FALSE)), times=mbt)
cat("Maximal eigensolution (negative matrix)\n")
print(agn)
cat("Geradin time=",msect(tgn$time),"sd=",msecr(tgn$time),"\n")


## ----timeger1, echo=TRUE-------------------------------------------------------------------------------------------------------------------------
naximp<-function(x, A=1){ # implicit moler A*x
   n<-length(x)
   y<-rep(0,n)
   for (i in 1:n){
      tt<-0.
      for (j in 1:n) {
          if (i == j) tt<-tt+i*x[i]
          else tt<-tt+(min(i,j) - 2)*x[j]
      }
      y[i]<- -tt # include negative sign
   }
   y
}

dyn.load("moler.so")
cat("Is the mat multiply loaded? ",is.loaded("moler"),"\n")

naxftn<-function(x, A) { # ignore second argument
   n<-length(x) # could speed up by having this passed
   vout<-rep(0,n) # purely for storage
   # NEED TO EXPLAIN -1 below
   res<-(-1)*(.Fortran("moler", n=as.integer(n), x=as.double(x), vout=as.double(vout)))$vout
}

require(microbenchmark)
nmax<-5
gtable<-matrix(NA, nrow=nmax, ncol=6) # to hold results
# loop over sizes
for (ni in 1:nmax){
  n<-100*ni
  x<-runif(n) # generate a vector 
  gtable[[ni, 1]]<-n
  AA<-molermat(n)
  BB<-diag(rep(1,n))
  tgax<-microbenchmark(ogax<-geradin(x, ax, bx, AA=-AA, BB=BB, control=list(trace=FALSE)), times=mbt)
  gtable[[ni, 2]]<-msect(tgax$time)
  tgaximp<-microbenchmark(ogaximp<-geradin(x, naximp, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)
  gtable[[ni, 3]]<-msect(tgaximp$time)
  tgaxftn<-microbenchmark(ogaxftn<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)
  gtable[[ni, 4]]<-msect(tgaxftn$time)
}

gtym<-data.frame(n=gtable[,1], ax=gtable[,2], aximp=gtable[,3], axftn=gtable[,4])
print(gtym)


# ## ----gerinr1, echo=TRUE, cache=TRUE--------------------------------------------------------------------------------------------------------------
# x<-runif(n)
# evalmax<-emax$evalmax
# evecmac<-emax$evecmax
# ogaxftn<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE))
# gvec<-ogaxftn$x
# gval<- -ogaxftn$RQ
# gvec<-sign(gvec[[1]])*gvec/sqrt(as.numeric(crossprod(gvec)))
# diff<-gvec-evecmax
# cat("Geradin eigenvalue - eigen result: ",gval-evalmax,"   max(abs(vector diff))=",
#        max(abs(diff)), "\n")



system("gfortran ./a25moler.f")
cat("Geradin fortran version a25moler.f")
tbld100<-msect(microbenchmark(AA<-molerfast(100), times=mbt)$time)
teig100<-msect(microbenchmark(a100<-eigen(AA), times=mbt)$time)
cat("eigen(): n=100 build time=",tbld100,"  eigen time=",teig100,"\n")
vecmin<-a100$vectors[,100]
vecmax<-a100$vectors[,1]
cat("eigen: Minimal Eigenvalue =",a100$values[100],"  RQ=", ray1(vecmin, AA),"\n")
tmin100<-msect(microbenchmark(system("./a.out <n100min.txt > out100min.txt"), times=mbt)$time)
res<-strsplit(readLines("out100min.txt"), " ")
resn<-as.numeric(res[[1]][which(res[[1]]!="")])
cat("A25RQM N=", as.integer(resn[1]), " matvec ops=", as.integer(resn[2]),
     " Min Est. EV=", resn[3], "  Gradient=",resn[4]," time=",msect(tmin100)," \n")
tmax100<-msect(microbenchmark(system("./a.out <n100max.txt > out100max.txt"), times=mbt)$time)
res<-strsplit(readLines("out100max.txt"), " ")
resn<-as.numeric(res[[1]][which(res[[1]]!="")])
cat("A25RQM N=", as.integer(resn[1]), " matvec ops=", as.integer(resn[2]),
     " Max Est. EV=", resn[3], "  Gradient=",resn[4]," time=",msect(tmax100)," \n")
cat("eigen: Maximal Eigenvalue =",a100$values[1],"  RQ=", ray1(vecmax, AA),"\n")


tbld200<-msect(microbenchmark(AA<-molerfast(200), times=mbt)$time)
teig200<-msect(microbenchmark(a200<-eigen(AA), times=mbt)$time)
cat("eigen(): n=200 build time=",tbld200,"  eigen time=",teig200,"\n")
vecmin<-a200$vectors[,200]
vecmax<-a200$vectors[,1]
cat("eigen: Minimal Eigenvalue =",a200$values[200],"  RQ=", ray1(vecmin, AA),"\n")
tmin200<-msect(microbenchmark(system("./a.out <n200min.txt > out200min.txt"), times=mbt)$time)
res<-strsplit(readLines("out200min.txt"), " ")
resn<-as.numeric(res[[1]][which(res[[1]]!="")])
cat("A25RQM N=", as.integer(resn[1]), " matvec ops=", as.integer(resn[2]),
     " Min Est. EV=", resn[3], "  Gradient=",resn[4]," time=",msect(tmin200)," \n")
tmax200<-msect(microbenchmark(system("./a.out <n200max.txt > out200max.txt"), times=mbt)$time)
res<-strsplit(readLines("out200max.txt"), " ")
resn<-as.numeric(res[[1]][which(res[[1]]!="")])
cat("A25RQM N=", as.integer(resn[1]), " matvec ops=", as.integer(resn[2]),
     " Max Est. EV=", resn[3], "  Gradient=",resn[4]," time=",msect(tmax200)," \n")
cat("eigen: Maximal Eigenvalue =",a200$values[1],"  RQ=", ray1(vecmax, AA),"\n")


tbld300<-msect(microbenchmark(AA<-molerfast(300), times=mbt)$time)
teig300<-msect(microbenchmark(a300<-eigen(AA), times=mbt)$time)
cat("eigen(): n=300 build time=",tbld300,"  eigen time=",teig300,"\n")
vecmin<-a300$vectors[,300]
vecmax<-a300$vectors[,1]
cat("eigen: Minimal Eigenvalue =",a300$values[300],"  RQ=", ray1(vecmin, AA),"\n")
tmin300<-msect(microbenchmark(system("./a.out <n300min.txt > out300min.txt"), times=mbt)$time)
res<-strsplit(readLines("out300min.txt"), " ")
resn<-as.numeric(res[[1]][which(res[[1]]!="")])
cat("A25RQM N=", as.integer(resn[1]), " matvec ops=", as.integer(resn[2]),
     " Min Est. EV=", resn[3], "  Gradient=",resn[4]," time=",msect(tmin300)," \n")
tmax300<-msect(microbenchmark(system("./a.out <n300max.txt > out300max.txt"), times=mbt)$time)
res<-strsplit(readLines("out300max.txt"), " ")
resn<-as.numeric(res[[1]][which(res[[1]]!="")])
cat("A25RQM N=", as.integer(resn[1]), " matvec ops=", as.integer(resn[2]),
     " Max Est. EV=", resn[3], "  Gradient=",resn[4]," time=",msect(tmax300)," \n")
cat("eigen: Maximal Eigenvalue =",a300$values[1],"  RQ=", ray1(vecmax, AA),"\n")


tbld400<-msect(microbenchmark(AA<-molerfast(400), times=mbt)$time)
teig400<-msect(microbenchmark(a400<-eigen(AA), times=mbt)$time)
cat("eigen(): n=400 build time=",tbld400,"  eigen time=",teig400,"\n")
vecmin<-a400$vectors[,400]
vecmax<-a400$vectors[,1]
cat("eigen: Minimal Eigenvalue =",a400$values[400],"  RQ=", ray1(vecmin, AA),"\n")
tmin400<-msect(microbenchmark(system("./a.out <n400min.txt > out400min.txt"), times=mbt)$time)
res<-strsplit(readLines("out400min.txt"), " ")
resn<-as.numeric(res[[1]][which(res[[1]]!="")])
cat("A25RQM N=", as.integer(resn[1]), " matvec ops=", as.integer(resn[2]),
     " Min Est. EV=", resn[3], "  Gradient=",resn[4]," time=",msect(tmin400)," \n")
tmax400<-msect(microbenchmark(system("./a.out <n400max.txt > out400max.txt"), times=mbt)$time)
res<-strsplit(readLines("out400max.txt"), " ")
resn<-as.numeric(res[[1]][which(res[[1]]!="")])
cat("A25RQM N=", as.integer(resn[1]), " matvec ops=", as.integer(resn[2]),
     " Max Est. EV=", resn[3], "  Gradient=",resn[4]," time=",msect(tmax400)," \n")
cat("eigen: Maximal Eigenvalue =",a400$values[1],"  RQ=", ray1(vecmax, AA),"\n")


tbld500<-msect(microbenchmark(AA<-molerfast(500), times=mbt)$time)
teig500<-msect(microbenchmark(a500<-eigen(AA), times=mbt)$time)
cat("eigen(): n=500 build time=",tbld500,"  eigen time=",teig500,"\n")
vecmin<-a500$vectors[,500]
vecmax<-a500$vectors[,1]
cat("eigen: Minimal Eigenvalue =",a500$values[500],"  RQ=", ray1(vecmin, AA),"\n")
tmin500<-msect(microbenchmark(system("./a.out <n500min.txt > out500min.txt"), times=mbt)$time)
res<-strsplit(readLines("out500min.txt"), " ")
resn<-as.numeric(res[[1]][which(res[[1]]!="")])
cat("A25RQM N=", as.integer(resn[1]), " matvec ops=", as.integer(resn[2]),
     " Min Est. EV=", resn[3], "  Gradient=",resn[4]," time=",msect(tmin500)," \n")
tmax500<-msect(microbenchmark(system("./a.out <n500max.txt > out500max.txt"), times=mbt)$time)
res<-strsplit(readLines("out500max.txt"), " ")
resn<-as.numeric(res[[1]][which(res[[1]]!="")])
cat("A25RQM N=", as.integer(resn[1]), " matvec ops=", as.integer(resn[2]),
     " Max Est. EV=", resn[3], "  Gradient=",resn[4]," time=",msect(tmax500)," \n")
cat("eigen: Maximal Eigenvalue =",a500$values[1],"  RQ=", ray1(vecmax, AA),"\n")



## ----persp1, echo=FALSE--------------------------------------------------------------------------------------------------------------------------
library(optimx)
cf<-data.frame(n=bmattym$n,
    eig=bmattym$eigentime,
    spg=spgtym$spgrqt,
    rcg=cgtime$tcgmin,
    ger=gtym$axftn )
eigen<-cf$eig/cf$ger
spg<-cf$spg/cf$ger
rcgmin<-cf$rcg/cf$ger
nsize<-cf$n
jn<-data.frame(nsize=nsize, eigen=eigen, spg=spg, rcgmin=rcgmin)
# joe<-write.csv(jn, file="jndata.csv")
plot(nsize,spg, pch=1, xlab="n", ylab="time ratio")
points(nsize, rcgmin, pch=3)
points(nsize, eigen, pch=4)
title("Ratio of eigensolution times to Geradin routine by matrix size")
points(nsize, rep(1,10), type="l")
#legend(50,70,c("spg", "rcgmin","eigen"), pch = c(1,3,4), lty = c(1,2,3))
legend(50,70,c("spg", "rcgmin","eigen"), pch = c(1,3,4))


## ----n2000a, echo=FALSE--------------------------------------------------------------------------------------------------------------------------
dyn.load("moler.so")
n<-2000
t2000b<-msect(microbenchmark(AA<-molerfast(n), times=mbt)$time)
t2000e<-msect(microbenchmark(evs<-eigen(AA), times=mbt)$time)
x<-runif(n)
t2000c<-msect(microbenchmark(ac<-optimr(x, fn=nobj, gr=ngrobj, method="ncg", 
                              AA=-AA), times=mbt)$time)
t2000g<-msect(microbenchmark(ag<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)$time)
cat("Times in seconds\n")
cat("Build =",t2000b," eigen():",t2000e,"  Rcgminu:", t2000c," Geradin:",t2000g,"\n")
cat("Ratios: build=", t2000b/t2000g, "eigen=",t2000e/t2000g,"  Rcgminu=",t2000c/t2000g,"\n")


## ----geradincodea, echo=TRUE---------------------------------------------------------------------------------------------------------------------
ax<-function(x, AA){
   u<-as.numeric(AA%*%x)
}
bx<-function(x, BB){
   v<-as.numeric(BB%*%x)
}
geradin<-function(x, ax, bx, AA, BB, control=list(trace=TRUE, maxit=1000)){
# Geradin minimize Rayleigh Quotient, Nash CMN Alg 25
# print(control)
  trace<-control$trace
  n<-length(x)
  tol<-n*n*.Machine$double.eps^2
  offset<-1e+5 # equality check offset
  if (trace) cat("geradin.R, using tol=",tol,"\n")
  ipr<-0 # counter for matrix mults
  pa<-.Machine$double.xmax
  R<-pa
  msg<-"no msg"
# step 1 -- main loop
  keepgoing<-TRUE
  while (keepgoing) {
    avec<-ax(x, AA); bvec<-bx(x, BB); ipr<-ipr+1
    xax<-as.numeric(crossprod(x, avec));  
    xbx<-as.numeric(crossprod(x, bvec));
    if (xbx <= tol) {
       keepgoing<-FALSE # not really needed
       msg<-"avoid division by 0 as xbx too small"
       break
    } 
    p0<-xax/xbx
    if (p0>pa) {
       keepgoing<-FALSE # not really needed
       msg<-"Rayleigh Quotient increased in step"
       break
    } 
    pa<-p0
    g<-2*(avec-p0*bvec)/xbx
    gg<-as.numeric(crossprod(g)) # step 6
    if (trace) cat("Before loop: RQ=",p0," after ",ipr," products, gg=",gg,"\n")
    if (gg<tol) { # step 7
       keepgoing<-FALSE # not really needed
       msg<-"Small gradient -- done"
       break
    } 
    t<- -g # step 8
    for (itn in 1:n) { # major loop step 9
       y<-ax(t, AA); z<-bx(t, BB); ipr<-ipr+1 # step 10
       tat<-as.numeric(crossprod(t, y)) # step 11
       xat<-as.numeric(crossprod(x, y)) 
       xbt<-as.numeric(crossprod(x, z)) 
       tbt<-as.numeric(crossprod(t, z)) 
       u<-tat*xbt-xat*tbt
       v<-tat*xbx-xax*tbt
       w<-xat*xbx-xax*xbt
       d<-v*v-4*u*w
       if (d<0) stop("Geradin: imaginary roots not possible") # step 13
       d<-sqrt(d) # step 14
       if (v>0) k<--2*w/(v+d) else k<-0.5*(d-v)/u
       xlast<-x # NOT as in CNM -- can be avoided with loop
       avec<-avec+k*y; bvec<-bvec+k*z # step 15, update
       x<-x+k*t
       xax<-xax+as.numeric(crossprod(x,avec))      
       xbx<-xbx+as.numeric(crossprod(x,bvec))      
       if (xbx<tol) stop("Geradin: xbx has become too small")
       chcount<-n - length(which((xlast+offset)==(x+offset)))
       if (trace) cat("Number of changed components = ",chcount,"\n")
       pn<-xax/xbx # step 17 different order
       if (chcount==0) {
         keepgoing<-FALSE # not really needed
         msg<-"Unchanged parameters -- done"
         break
       }
       if (pn >= p0) {
         if (trace) cat("RQ not reduced, restart\n")
         break # out of itn loop, not while loop (TEST!)
       }
       p0<-pn # step 19
       g<-2*(avec-pn*bvec)/xbx
       gg<-as.numeric(crossprod(g))
       if (trace) cat("Itn", itn," RQ=",p0," after ",ipr," products, gg=",gg,"\n")
       if (gg<tol){ # step 20
         if (trace) cat("Small gradient in iteration, restart\n")
         break # out of itn loop, not while loop (TEST!)
       }
       xbt<-as.numeric(crossprod(x,z)) # step 21
       w<-y-pn*z # step 22
       tabt<-as.numeric(crossprod(t,w))
       beta<-as.numeric(crossprod(g,(w-xbt*g)))
       beta<-beta/tabt # step 23
       t<-beta*t-g
    } # end loop on itn -- step 24
  } # end main loop -- step 25
  ans<-list(x=x, RQ=p0, ipr=ipr, msg=msg) # step 26
}

