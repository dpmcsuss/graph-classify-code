
#  see <http://mit.edu/lrv/www/elegans/>
#  for labeled graphs of c. elegans (in matlab).
#  ---cep

# We consider k(n) nrst nbr classification
# of labeled multigraphs on 279 vertices,
# under Frobenius norm.
# Class-conditional distributions F_{G|Y=j}
# are specified via Poisson parameters for each edge.
# Class-conditional distribution F_{G|Y=0}
# is given by "celegans + noise".
# (noise parameter "eps".)
# Class-conditional distribution F_{G|Y=1}
# is given by "celegans + noise + egg".
# (noise parameter "eps", egg parameters "m" & "sig".)
# NB: a minor but deterministic class-conditional difference
# is imposed, so that L^{\star} = 0.

getLhat = function()
{
load("celegansGraph.Rd")  ## NB: dim(Ac) => 279 279
set.seed(12345)
tstn = 500  ## NB: L=0.1 & 900 test observations => stdev(Lhat^(R)) = 0.01
maxn = 200  ## maximum number of training observations per class
## NB: we've conditioned on cc ss n_0=n_1 for both training and testing
n = 279
m = 14  ## egg consists of 5% of the vertices ... just because.
thesem = sample(n,m)  ## NB: egg vertices; same for all mc replicates!
eps = 0.05
sig = 5
E0 = Ac + matrix(runif(n^2,0,eps),nrow=n,ncol=n)
E1=E0
egg = matrix(runif(m^2,-sig/2,sig/2),nrow=m,ncol=m)
E1[thesem,thesem] = E1[thesem,thesem] + egg
for(i in thesem) for(j in thesem) if(E1[i,j]<0) E1[i,j]=0
## NB: Frobenius norm => we can use vector representation for graphs!
Amat = matrix(0,nrow=(2*maxn),ncol=(n*n))  ## training data
Bmat = matrix(0,nrow=(2*tstn),ncol=(n*n))  ## test data
for(i in 1:n) for(j in 1:n)
 {
 indexij = i+(j-1)*n
 Amat[1:maxn,indexij]            = rpois(maxn,E0[i,j]) ## train, class 0
 Amat[(maxn+1):(2*maxn),indexij] = rpois(maxn,E1[i,j]) ## train, class 1
 Bmat[1:tstn,indexij]            = rpois(tstn,E0[i,j]) ## test, class 0
 Bmat[(tstn+1):(2*tstn),indexij] = rpois(tstn,E1[i,j]) ## test, class 1
 }
# NB: a minor but deterministic class-conditional difference is imposed, so that L^*=0.
Amat[1:maxn,8495]=1 ; Amat[(maxn+1):(2*maxn),8495]=2 ## trick so L^*=0
Bmat[1:tstn,8495]=1 ; Bmat[(tstn+1):(2*tstn),8495]=2 ## trick so L^*=0
trueclasslabel=c(rep(0,tstn),rep(1,tstn))
dm = matrix(0,nrow=(2*tstn),ncol=(2*maxn))
for(i in 1:(2*tstn))
 dm[i,] = apply(Amat,1,function(x) sqrt(sum((Bmat[i,]-x)^2)))
oi = matrix(0,nrow=(2*tstn),ncol=(2*maxn))
for(i in 1:(2*tstn))
 oi[i,] = order(dm[i,])

k=Lhat=NULL
nvec=seq(4,maxn,by=4)
nvec=seq(5,150,by=5)
for(z in 1:length(nvec))
 {
 trainn = nvec[z]
 trainobs = c( 1:trainn , (maxn+1):(maxn+trainn) )
 #k[z]=floor(log(trainn)) ; if( (k[z]%%2)==0 ) k[z]=k[z]+1
 k[z]=floor(sqrt(trainn)) ; if( (k[z]%%2)==0 ) k[z]=k[z]-1
 knnclasslabel=rep(1,(2*tstn))
 for(i in 1:(2*tstn))
  if(sum( oi[i,trainobs][1:k[z]] <= trainn )>(k[z]/2)) knnclasslabel[i]=0
 Lhat[z] = sum(trueclasslabel != knnclasslabel)/(2*tstn)
 print(paste(nvec[z],k[z],round(Lhat[z],3)))
 }

# "c elegans graph classification simulation results"
plot(nvec,Lhat,pch=19,col=(k-1)/2,ylim=c(0,1),type="b",xlab=expression(n[j]),ylab=expression(hat(L)))
# colors are (k(n)-1)/2, where k(n)=sqrt(n/2) => uc.
# standard errors are sqrt(L*(1-L)/1000).
# Eg, n_j = 150 ; k(n) = 11 ; Lhat = 0.069.
# standard errors are less than 0.01 at n_j = 150;
# we reject H0: L^{\star} >= 0.10 at alpha=0.05.
# NB: L^{\star} = 0.

par(pty="s")
pdf("Lhatplot.pdf")
plot(nvec,Lhat,pch=19,col=(k-1)/2,ylim=c(0,1),type="b",xlab=expression(n[j]),ylab=expression(hat(L)))
dev.off()
