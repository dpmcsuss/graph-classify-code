
# see <http://mit.edu/lrv/www/elegans/> for labeled graphs of c. elegans.

# We consider classification
# of labeled multigraphs (directed, with loops)
# on 279 vertices.
# Class-conditional distributions F_{G|Y=j}
# are specified via Poisson parameters for each edge.
# Class-conditional distribution F_{G|Y=0}
# is given by "celegans + noise".
# (noise parameter "eps".)
# Class-conditional distribution F_{G|Y=1}
# is given by "celegans + noise + egg".
# (noise parameter "eps", egg parameter "sig".)
# egg vertices:
#  based on Chalasni et al.,
#  the following neurons are involved in odor behaviors:
#  AWC, AIA, AIB, AIZ, AIY.
#  each of the above actually has a left neuron and a right neuron,
#  meaning 10 neurons total.
#  so, in our sim, the egg is defined by this set of m=10 neurons.

#load("~/Research/necog/supervenience/data/celegansGraph.Rd")
#A0=A0[1:71,1:71]
#n=sqrt(length(A0))

rm(list=ls())
A0=read.delim(file("~/Research/necog/BLSA1/data/john/MeanFA.csv"), header=FALSE,sep="")
n = length(A0)
A0=as.matrix(A0,nrow=n,ncol=n)
A0=matrix(A0,nrow=n,ncol=n)
A0=A0+t(A0)
NaNind=is.na(A0)
A0[NaNind]=0
A0=A0*10




set.seed(12345)
tstn = 500  ## NB: L=0.1 & 1000 test observations => stdev(Lhat^(R)) < 0.01
maxn = 200  ## maximum number of training observations per class
## NB: we condition on class-cond'l sample size n_0=n_1 for both training and testing
##  m=10 and
##  AWCL=69 AIBL=80 AWCR=82 AIRB=94 AIAL=110 AIAR=127 AIZR=129 AIZL=133 AIYR=134 AIYL=138
# m=10
# thesem = c(69,80,82,94,110,127,129,133,134,138)
m=5
thesem=1:5
eps = 0.05  # noise parameter "eps"
sig = 4     # egg parameter "sig"
E0 = A0 + matrix(runif(n^2,0,eps),nrow=n,ncol=n)
E1=E0
egg = matrix(runif(m^2,-sig/2,sig/2),nrow=m,ncol=m)
E1[thesem,thesem] = E1[thesem,thesem] + egg
for(i in thesem) for(j in thesem) if(E1[i,j]<0) E1[i,j]=0

Amatlist = Bmatlist = NULL
for(i in 1:maxn)
 Amatlist[[i]] = matrix(rpois(n^2,E0),nrow=n,ncol=n)
for(i in (maxn+1):(2*maxn))
 Amatlist[[i]] = matrix(rpois(n^2,E1),nrow=n,ncol=n)
for(i in 1:tstn)
 Bmatlist[[i]] = matrix(rpois(n^2,E0),nrow=n,ncol=n)
for(i in (tstn+1):(2*tstn))
 Bmatlist[[i]] = matrix(rpois(n^2,E1),nrow=n,ncol=n)

nvec=seq(5,180,by=5)


LhatNB=Lhati=Lhatc=NULL
for(z in 1:length(nvec))
 {
 trainn = nvec[z]
beta=1
lambdahat0 = matrix(0,nrow=n,ncol=n)
for(i in 1:trainn)
 lambdahat0 = lambdahat0 + Amatlist[[i]]
lambdahat0 = lambdahat0 / trainn
lambdahat1 = matrix(0,nrow=n,ncol=n)
for(i in (maxn+1):(maxn+trainn))
 lambdahat1 = lambdahat1 + Amatlist[[i]]
lambdahat1 = lambdahat1 / trainn
 for(i in 1:length(lambdahat0)) if(lambdahat0[i]==0)
lambdahat0[i]=1/(2*trainn)
 for(i in 1:length(lambdahat1)) if(lambdahat1[i]==0)
lambdahat1[i]=1/(2*trainn)
deltahat = abs(lambdahat1-lambdahat0)

# par(mfrow=c(3,1))

## assume incoherent
#plot(c(deltahat))
#points(sigdims,c(deltahat[thesem,thesem]),col=2,pch=19)
# plot(rev(sort(deltahat))[1:100])
mi = matrix(1:n^2,ncol=n,nrow=n,byrow=F)
sigdims = c(mi[thesem,thesem])
wwwi = which(order(deltahat,decreasing=T) %in% sigdims)
# points( wwwi,rev(sort(deltahat))[wwwi],col=2,pch=19)

## assume semicoherent
deg = apply(deltahat,1,sum) + apply(deltahat,2,sum)
#plot(deg)
#points(thesem,deg[thesem],col=2,pch=19)
# plot(rev(sort(deg))[1:100])
wwws = which(order(deg,decreasing=T) %in% thesem)
# points( wwws,rev(sort(deg))[wwws],col=2,pch=19)

## assume coherent
deltahateps = deltahat*(deltahat>beta)
loc=NULL
for(i in 1:n)
 {
 nbhd = c(i,which(deltahateps[i,]>0))
 loc1 = sum(deltahateps[nbhd,nbhd])
 nbhd = c(i,which(deltahateps[,i]>0))
 loc2 = sum(deltahateps[nbhd,nbhd])
 loc[i]=loc1+loc2
 }
#plot(loc)
#points(thesem,loc[thesem],col=2,pch=19)
# plot(rev(sort(loc))[1:100])
wwwc = which(order(loc,decreasing=T) %in% thesem)
# points( wwwc,rev(sort(loc))[wwwc],col=2,pch=19)

# dev.copy(pdf,"subspacerecovery.pdf")
# dev.off()

# classify, using all n^2==77841 dimensions:
 trueclasslabel=c(rep(0,tstn),rep(1,tstn))
 iclasslabel=rep(1,(2*tstn))
 for(i in 1:(2*tstn))
 if(sum(log(dpois(c(Bmatlist[[i]]),c(lambdahat0)))) >
sum(log(dpois(c(Bmatlist[[i]]),c(lambdahat1))))) iclasslabel[i]=0
 LhatNB[z] = sum(trueclasslabel != iclasslabel)/(2*tstn)
print(paste("Lhat-NaiveBayes",LhatNB[z]))

# classify, using 100 dimensions recovered via incoherent assumption:
sigdimesti = order(deltahat,decreasing=T)[1:100]
 trueclasslabel=c(rep(0,tstn),rep(1,tstn))
 iclasslabel=rep(1,(2*tstn))
 for(i in 1:(2*tstn))
 if(sum(log(dpois(c(Bmatlist[[i]])[sigdimesti],c(lambdahat0)[sigdimesti])))
>
sum(log(dpois(c(Bmatlist[[i]])[sigdimesti],c(lambdahat1)[sigdimesti]))))
iclasslabel[i]=0
 Lhati[z] = sum(trueclasslabel != iclasslabel)/(2*tstn)
print(paste("Lhat-incoherentsubspace",Lhati[z]))

# classify, using 100 dimensions recovered via coherent assumption:
vest = order(loc,decreasing=T)[1:10]
sigdimest = c(mi[vest,vest])
 trueclasslabel=c(rep(0,tstn),rep(1,tstn))
 classlabel=rep(1,(2*tstn))
 for(i in 1:(2*tstn))
 if(sum(log(dpois(c(Bmatlist[[i]])[sigdimest],c(lambdahat0)[sigdimest])))
> sum(log(dpois(c(Bmatlist[[i]])[sigdimest],c(lambdahat1)[sigdimest]))))
classlabel[i]=0
 Lhatc[z] = sum(trueclasslabel != classlabel)/(2*tstn)
print(paste("Lhat-coherentsubspace",Lhatc[z]))
}

# NB: McNemar's test => statistical significance for
# NB->incoherent->coherent improvement


# kNN
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
trueclasslabel=c(rep(0,tstn),rep(1,tstn))
dm = matrix(0,nrow=(2*tstn),ncol=(2*maxn))
for(i in 1:(2*tstn))
 dm[i,] = apply(Amat,1,function(x) sqrt(sum((Bmat[i,]-x)^2)))
oi = matrix(0,nrow=(2*tstn),ncol=(2*maxn))
for(i in 1:(2*tstn))
 oi[i,] = order(dm[i,])


k=Lhat=NULL
nvec=seq(5,180,by=5)
for(z in 1:length(nvec))
 {
 trainn = nvec[z]
 trainobs = c( 1:trainn , (maxn+1):(maxn+trainn) )
 k[z]=floor(sqrt(16*trainn)) ; if( (k[z]%%2)==0 ) k[z]=k[z]-1
 knnclasslabel=rep(1,(2*tstn))
 for(i in 1:(2*tstn))
  if(sum( oi[i,trainobs][1:k[z]] <= trainn )>(k[z]/2)) knnclasslabel[i]=0
 Lhat[z] = sum(trueclasslabel != knnclasslabel)/(2*tstn)
print(paste("Lhat-knn",Lhat[z]))
 }














