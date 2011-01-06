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

load("celegansGraph.Rd")  ## NB: dim(Ac) => 279 279
Ac=Ac[1:70,1:70];
set.seed(12345)
tstn = 500  ## NB: L=0.1 & 1000 test observations => stdev(Lhat^(R)) < 0.01
maxn = 200  ## maximum number of training observations per class
## NB: we condition on class-cond'l sample size n_0=n_1 for both training and testing
n = 70
##  m=10 and
##  AWCL=69 AIBL=80 AWCR=82 AIRB=94 AIAL=110 AIAR=127 AIZR=129 AIZL=133 AIYR=134 AIYL=138
m=5
thesem=1:m
eps = 0.05  # noise parameter "eps"
sig = 5     # egg parameter "sig"
E0 = Ac + matrix(runif(n^2,0,eps),nrow=n,ncol=n)
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


trainn=5
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

par(mfrow=c(3,1))

## assume incoherent
#plot(c(deltahat))
#points(sigdims,c(deltahat[thesem,thesem]),col=2,pch=19)
plot(rev(sort(deltahat))[1:100])
mi = matrix(1:n^2,ncol=n,nrow=n,byrow=F)
sigdims = c(mi[thesem,thesem])
wwwi = which(order(deltahat,decreasing=T) %in% sigdims)
points( wwwi,rev(sort(deltahat))[wwwi],col=2,pch=19)

## assume semicoherent
deg = apply(deltahat,1,sum) + apply(deltahat,2,sum)
#plot(deg)
#points(thesem,deg[thesem],col=2,pch=19)
plot(rev(sort(deg))[1:100])
wwws = which(order(deg,decreasing=T) %in% thesem)
points( wwws,rev(sort(deg))[wwws],col=2,pch=19)

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
plot(rev(sort(loc))[1:100])
wwwc = which(order(loc,decreasing=T) %in% thesem)
points( wwwc,rev(sort(loc))[wwwc],col=2,pch=19)

dev.copy(pdf,"subspacerecovery.pdf")
dev.off()
# > getwd()
# [1] "C:/Users/priebe/Documents/My Dropbox/Desktop0729/Projects/Vogelstein/PhilosophyPaper/c elegans"
# workingDirectory <-"/Users/priebe/Documents/My Dropbox/Desktop0729/Projects/sss-Andrey Rukhin"
# setwd(workingDirectory)
# dev.copy(pdf,paste(outputDirectory,"Size_m50_q0pt5.pdf",sep=""))
# dev.off()

# classify, using all n^2==77841 dimensions:
 trueclasslabel=c(rep(0,tstn),rep(1,tstn))
 iclasslabel=rep(1,(2*tstn))
 for(i in 1:(2*tstn))
 if(sum(log(dpois(c(Bmatlist[[i]]),c(lambdahat0)))) >
sum(log(dpois(c(Bmatlist[[i]]),c(lambdahat1))))) iclasslabel[i]=0
 LhatNB = sum(trueclasslabel != iclasslabel)/(2*tstn)
print(paste("Lhat-NaiveBayes",LhatNB))
# [1] "Lhat-NaiveBayes 0.039"

# classify, using 100 dimensions recovered via incoherent assumption:
sigdimesti = order(deltahat,decreasing=T)[1:100]
 # > sum(wwwi<=100)
 # [1] 16
 #  => 16 signal dimensions recovered
 trueclasslabel=c(rep(0,tstn),rep(1,tstn))
 iclasslabel=rep(1,(2*tstn))
 for(i in 1:(2*tstn))
 if(sum(log(dpois(c(Bmatlist[[i]])[sigdimesti],c(lambdahat0)[sigdimesti])))

sum(log(dpois(c(Bmatlist[[i]])[sigdimesti],c(lambdahat1)[sigdimesti]))))
iclasslabel[i]=0
 Lhati = sum(trueclasslabel != iclasslabel)/(2*tstn)
print(paste("Lhat-incoherentsubspace",Lhati))
# [1] "Lhat-incoherentsubspace 0.013"

# classify, using 100 dimensions recovered via coherent assumption:
vest = order(loc,decreasing=T)[1:10]
sigdimest = c(mi[vest,vest])
 # > sum(wwwc<=10)
 # [1] 8
 #  => 64 signal dimensions recovered
 trueclasslabel=c(rep(0,tstn),rep(1,tstn))
 classlabel=rep(1,(2*tstn))
 for(i in 1:(2*tstn))
 if(sum(log(dpois(c(Bmatlist[[i]])[sigdimest],c(lambdahat0)[sigdimest])))
 sum(log(dpois(c(Bmatlist[[i]])[sigdimest],c(lambdahat1)[sigdimest]))))
classlabel[i]=0
 Lhatc = sum(trueclasslabel != classlabel)/(2*tstn)
print(paste("Lhat-coherentsubspace",Lhatc))
# [1] "Lhat-coherentsubspace 0"

# NB: McNemar's test => statistical significance for
# NB->incoherent->coherent improvement