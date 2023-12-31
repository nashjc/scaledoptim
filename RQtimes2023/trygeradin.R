rq<-  as.numeric(crossprod(x, crossprod(AA,x)))
}
# Chunk 12: raynum3
ray3<-function(x, AA, ax=axftn){
# ax is a function to form AA%*%x
rq<- - as.numeric(crossprod(x, ax(x, AA)))
}
# Chunk 13: axm
ax<-function(x, AA){
u<- as.numeric(AA%*%x)
}
axx<-function(x, AA){
u<- as.numeric(crossprod(AA, x))
}
# Chunk 14: aximp
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
# Chunk 15: axmfcode
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
# Chunk 16: axftn
dyn.load("moler.so")
cat("Is the mat multiply loaded? ",is.loaded("moler"),"\n")
axftn<-function(x, AA=1) { # ignore second argument
n<-length(x) # could speed up by having this passed
vout<-rep(0,n) # purely for storage
res<-(.Fortran("moler", n=as.integer(n), x=as.double(x), vout=as.double(vout)))$vout
}
# Chunk 17: timeax1
dyn.load("moler.so")
cat("Is the mat multiply loaded? ",is.loaded("moler"),"\n")
require(microbenchmark)
nmax<-10
ptable<-matrix(NA, nrow=nmax, ncol=6) # to hold results
# loop over sizes
for (ni in 1:nmax){
n<-50*ni
x<-runif(n) # generate a vector
ptable[[ni, 1]]<-n
AA<-molermat(n)
tax<-system.time(oax<-replicate(100,ax(x, AA))[,1])[[1]]
taxx<-system.time(oaxx<-replicate(100,axx(x, AA))[,1])[[1]]
if (! identical(oax, oaxx)) stop("oaxx NOT correct")
taxftn<-system.time(oaxftn<-replicate(100,axftn(x, AA=1))[,1])[[1]]
if (! identical(oax, oaxftn)) stop("oaxftn NOT correct")
taximp<-system.time(oaximp<-replicate(100,aximp(x, AA=1))[,1])[[1]]
if (! identical(oax, oaximp)) stop("oaximp NOT correct")
taxmfi<-system.time(oaxmfi<-replicate(100,axmolerfast(x, AA=1))[,1])[[1]]
if (! identical(oax, oaxmfi)) stop("oaxmfi NOT correct")
ptable[[ni, 2]]<-tax
ptable[[ni, 3]]<-taxx
ptable[[ni, 4]]<-taxftn
ptable[[ni, 5]]<-taximp
ptable[[ni, 6]]<-taxmfi
}
axtym<-data.frame(n=ptable[,1], ax=ptable[,2], axx=ptable[,3],  axftn=ptable[,4],
aximp=ptable[,5], axmfast=ptable[,6])
print(axtym)
# Chunk 18: adjaxtime
adjtym<-data.frame(n=axtym$n, axx1=axtym$axx+1*bmattym$buildi,
axxz=axtym$axx+100*bmattym$buildi,
axftn=axtym$axftn, aximp=axtym$aximp)
print(adjtym)
# Chunk 19: rqtime1
dyn.load("moler.so")
n<-500
x<-runif(n) # generate a vector
AA<-molermat(n)
tdi<-microbenchmark(rdi<-rqdir(x, AA))$time
cat("Direct algorithm: ",mean(tdi)*0.001,"\n")
t1i<-microbenchmark(r1i<-ray1(x, AA))$time
cat("ray1: mat-mult algorithm: ", mean(t1i)*0.001,"\n")
t2i<-microbenchmark(r2i<-ray2(x, AA))$time
cat("ray2: crossprod algorithm: ",mean(t2i)*0.001,"\n")
t3fi<-microbenchmark(r3i<-ray3(x, AA, ax=axftn))$time
cat("ray3: ax Fortran + crossprod: ",mean(t3fi)*0.001,"\n")
t3ri<-microbenchmark(r3i<-ray3(x, AA, ax=axmolerfast))$time
cat("ray3: ax fast R implicit + crossprod: ",mean(t3ri)*0.001,"\n")
# Chunk 20: rayspg1
rqt<-function(x, AA){
rq<-as.numeric(crossprod(x, crossprod(AA,x)))
}
proj<-function(x) { sign(x[1])*x/sqrt(crossprod(x)) }
require(BB)
n<-100
x<-rep(1,n)
AA<-molermat(n)
evs<-eigen(AA)
tmin<-microbenchmark(amin<-spg(x, fn=rqt, project=proj, control=list(trace=FALSE), AA=AA), times=mbt)$time
#amin
tmax<-microbenchmark(amax<-spg(x, fn=rqt, project=proj, control=list(trace=FALSE), AA=-AA), times=mbt)$time
#amax
evalmax<-evs$values[1]
evecmax<-evs$vectors[,1]
evecmax<-sign(evecmax[1])*evecmax/sqrt(as.numeric(crossprod(evecmax)))
emax<-list(evalmax=evalmax, evecmax=evecmax)
# save(emax, file="temax.Rdata")
evalmin<-evs$values[n]
evecmin<-evs$vectors[,n]
evecmin<-sign(evecmin[1])*evecmin/sqrt(as.numeric(crossprod(evecmin)))
avecmax<-amax$par
avecmin<-amin$par
avecmax<-sign(avecmax[1])*avecmax/sqrt(as.numeric(crossprod(avecmax)))
avecmin<-sign(avecmin[1])*avecmin/sqrt(as.numeric(crossprod(avecmin)))
cat("minimal eigensolution: Value=",amin$value,"in time ",tmin,"\n")
cat("Eigenvalue - result from eigen=",amin$value-evalmin,"  vector max(abs(diff))=",
max(abs(avecmin-evecmin)),"\n\n")
#print(amin$par)
cat("maximal eigensolution: Value=",-amax$value,"in time ",tmax,"\n")
cat("Eigenvalue - result from eigen=",-amax$value-evalmax,"  vector max(abs(diff))=",
max(abs(avecmax-evecmax)),"\n\n")
#print(amax$par)
# Chunk 21: runspg2
require(compiler)
require(BB)
nmax<-10
stable<-matrix(NA, nrow=nmax, ncol=3) # to hold results
# loop over sizes
for (ni in 1:nmax){
n<-50*ni
x<-runif(n) # generate a vector
AA<-molermat(n) # make sure defined
stable[[ni, 1]]<-n
tbld<-microbenchmark(AA<-molerfast(n), times=mbt)
tspg<-microbenchmark(aspg<-spg(x, fn=rqt, project=proj, control=list(trace=FALSE), AA=-AA), times=mbt)
stable[[ni, 2]]<-mean(tspg$time)*0.001
stable[[ni, 3]]<-mean(tbld$time)*0.001
}
spgtym<-data.frame(n=stable[,1], spgrqt=stable[,2], tbld=stable[,3])
print(round(spgtym,0))
# Chunk 22: runopx1
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
# require(optplus)
# mset<-c("L-BFGS-B", "BFGS", "CG", "spg", "ucminf", "nlm", "nlminb", "Rvmmin", "Rcgmin")
mset<-c("L-BFGS-B", "BFGS", "ncg", "spg", "ucminf", "nlm", "nlminb", "nvm")
nmax<-5
for (ni in 1:nmax){
n<-20*ni
x<-runif(n) # generate a vector
AA<-molerfast(n) # make sure defined
aall<-opm(x, fn=nobj, gr=ngrobj, method=mset, AA=-AA,
control=list(starttests=FALSE, dowarn=FALSE))
# optansout(aall, NULL)
summary(aall, order=value, )
cat("Above for n=",n," \n")
}
# Chunk 23: rcgrun1
ctable<-matrix(NA, nrow=10, ncol=2)
nmax<-10
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
# Chunk 24: geradincode
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
# Chunk 25: rungeradin1
cat("Test geradin with explicit matrix multiplication\n")
n<-10
AA<-molermat(n)
BB=diag(rep(1,n))
x<-runif(n)
tg<-microbenchmark(ag<-geradin(x, ax, bx, AA=AA, BB=BB,
control=list(trace=FALSE)), times=mbt)
cat("Minimal eigensolution\n")
print(ag)
cat("Geradin time=",mean(tg$time),"\n")
tgn<-microbenchmark(agn<-geradin(x, ax, bx, AA=-AA, BB=BB,
control=list(trace=FALSE)), times=mbt)
cat("Maximal eigensolution (negative matrix)\n")
print(agn)
cat("Geradin time=",mean(tgn$time),"\n")
# Chunk 26: timeger1
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
res<-(-1)*(.Fortran("moler", n=as.integer(n), x=as.double(x), vout=as.double(vout)))$vout
}
require(microbenchmark)
nmax<-10
gtable<-matrix(NA, nrow=nmax, ncol=6) # to hold results
# loop over sizes
for (ni in 1:nmax){
n<-50*ni
x<-runif(n) # generate a vector
gtable[[ni, 1]]<-n
AA<-molermat(n)
BB<-diag(rep(1,n))
tgax<-microbenchmark(ogax<-geradin(x, ax, bx, AA=-AA, BB=BB, control=list(trace=FALSE)), times=mbt)
gtable[[ni, 2]]<-mean(tgax$time)
tgaximp<-microbenchmark(ogaximp<-geradin(x, naximp, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)
gtable[[ni, 3]]<-mean(tgaximp$time)
tgaxftn<-microbenchmark(ogaxftn<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)
gtable[[ni, 4]]<-mean(tgaxftn$time)
}
gtym<-data.frame(n=gtable[,1], ax=gtable[,2], aximp=gtable[,3], axftn=gtable[,4])
print(gtym)
# Chunk 27: gercheck1
n<-100
x<-runif(n)
# emax<-load("temax.Rdata")
evalmax<-emax$evalmax
evecmac<-emax$evecmax
ogaxftn<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE))
gvec<-ogaxftn$x
gval<- -ogaxftn$RQ
gvec<-sign(gvec[[1]])*gvec/sqrt(as.numeric(crossprod(gvec)))
diff<-gvec-evecmax
cat("Geradin diff eigenval from eigen result: ",gval-evalmax,"   max(abs(vector diff))=",
max(abs(diff)), "\n")
cgtime
eigen<-cf$eig/cf$ger
cf<-data.frame(n=bmattym$n,
eig=bmattym$eigentime,
spg=spgtym$spgrqt,
rcg=cgtime$tcgmin,
ger=gtym$axftnc )
gtym
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
joe<-write.csv(jn, file="jndata.csv")
plot(nsize,spg, pch=1, xlab="n", ylab="time ratio")
points(nsize, rcgmin, pch=3)
points(nsize, eigen, pch=4)
title("Ratio of eigensolution times to Geradin routine by matrix size")
points(nsize, rep(1,10), type="l")
#legend(50,70,c("spg", "rcgmin","eigen"), pch = c(1,3,4), lty = c(1,2,3))
legend(50,70,c("spg", "rcgmin","eigen"), pch = c(1,3,4))
t2000g<-microbenchmark(ag<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)
dyn.load("moler.so")
n<-2000
t2000b<-microbenchmark(AA<-molerfast(n), times=mbt)
t2000e<-microbenchmark(evs<-eigen(AA), times=mbt)
x<-runif(n)
t2000c<-microbenchmark(ac<-Rcgminu(x, fn=nobj, gr=ngrobj, AA=-AA), times=mbt)
t2000g<-microbenchmark(ag<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)
cat("Times in seconds\n")
cat("Build =",t2000b," eigen():",t2000e,"  Rcgminu:", t2000c," Geradin:",t2000g,"\n")
dyn.load("moler.so")
n<-2000
t2000b<-microbenchmark(AA<-molerfast(n), times=mbt)
t2000e<-microbenchmark(evs<-eigen(AA), times=mbt)
x<-runif(n)
t2000c<-microbenchmark(ac<-Rcgminu(x, fn=nobj, gr=ngrobj, AA=-AA), times=mbt)
t2000g<-mean(microbenchmark(ag<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)$time)
cat("Times in seconds\n")
cat("Build =",t2000b," eigen():",t2000e,"  Rcgminu:", t2000c," Geradin:",t2000g,"\n")
dyn.load("moler.so")
n<-2000
t2000b<-microbenchmark(AA<-molerfast(n), times=mbt)
t2000e<-microbenchmark(evs<-eigen(AA), times=mbt)
x<-runif(n)
t2000c<-mean(microbenchmark(ac<-Rcgminu(x, fn=nobj, gr=ngrobj, AA=-AA), times=mbt)$time)
t2000g<-mean(microbenchmark(ag<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)$time)
cat("Times in seconds\n")
cat("Build =",t2000b," eigen():",t2000e,"  Rcgminu:", t2000c," Geradin:",t2000g,"\n")
dyn.load("moler.so")
n<-2000
t2000b<-microbenchmark(AA<-molerfast(n), times=mbt)
t2000e<-microbenchmark(evs<-eigen(AA), times=mbt)
x<-runif(n)
t2000c<-mean(microbenchmark(ac<-optimr(x, fn=nobj, gr=ngrobj, method="ncg",
AA=-AA), times=mbt)$time)
t2000g<-mean(microbenchmark(ag<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)$time)
cat("Times in seconds\n")
cat("Build =",t2000b," eigen():",t2000e,"  Rcgminu:", t2000c," Geradin:",t2000g,"\n")
dyn.load("moler.so")
n<-2000
t2000b<-mean(microbenchmark(AA<-molerfast(n), times=mbt)$time)
t2000e<-mean(microbenchmark(evs<-eigen(AA), times=mbt)$time)
x<-runif(n)
t2000c<-mean(microbenchmark(ac<-optimr(x, fn=nobj, gr=ngrobj, method="ncg",
AA=-AA), times=mbt)$time)
t2000g<-mean(microbenchmark(ag<-geradin(x, naxftn, ident, AA=1, BB=1, control=list(trace=FALSE)), times=mbt)$time)
cat("Times in seconds\n")
cat("Build =",t2000b," eigen():",t2000e,"  Rcgminu:", t2000c," Geradin:",t2000g,"\n")
cat("Ratios: build=", t2000b/t2000g, "eigen=",t2000e/t2000g,"  Rcgminu=",t2000c/t2000g,"\n")
dir()
library(microbenchmark)
system.time(system("a.out < n500.txt"))
microbenchmark(system("a.out < n500.txt"))
microbenchmark(system("./a.out < n500.txt"))
microbenchmark(system("./a.out < n500.txt > tmp.txt"))
mbt<-25
ta25<-microbenchmark(system("./a.out < n500.txt > tmp.txt"), times=mbt)
mean(ta256$times)
mean(ta25$times)
ta25
mean(ta25$time)
mean(ta25$time)*.001
mbt<-25 # default to 5 repetitions in microbenchmark while sorting out text
msect<-function(times){
round(mean(times)*0.001,0)
}
msecr<-function(times){
#   round((max(times)-min(times))*0.001,0)
round(sd(times)*0.001,0)
}
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
require(microbenchmark)
nmax<-10
mtable<-matrix(NA, nrow=nmax, ncol=5) # to hold results
rtable<-mtable
# loop over sizes
for (ni in 1:nmax){
n<-50*ni
mtable[[ni, 1]]<-n
rtable[[ni, 1]]<-n
# Note "unit" argument is ONLY for display. time is in nanoseconds
ti<-microbenchmark(ai<-molermat(n), unit="us", times=mbt)$time
tfi<-microbenchmark(afi<-molerfast(n), unit="us", times=mbt)$time
if (! identical(ai, afi)) stop("Different outcomes == molermat, molerfast")
osize<-object.size(ai)
tevs<-microbenchmark(evs<-eigen(ai), unit="us", times=mbt)$time
mtable[[ni,2]]<-msect(ti)
mtable[[ni,3]]<-osize
mtable[[ni,4]]<-msect(tevs)
mtable[[ni,5]]<-msect(tfi)
rtable[[ni,2]]<-msecr(ti)
rtable[[ni,3]]<-osize
rtable[[ni,4]]<-msecr(tevs)
rtable[[ni,5]]<-msecr(tfi)
}
AA<-molerfast(500)
system.time(eigen(AA))
.080*1e6
ls()
ta25
teig<-microbenchmark(eigen(AA), times=mbt)
msect(teig$time)
msect(microbenchmark(system("./a.out < n500.txt >tmp.txt"),times=mbt)$time)
savehistory("trygeradin")
