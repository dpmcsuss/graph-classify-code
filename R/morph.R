 # n0,n1,n,d are global variables
 n0=20
 n1=150
 n=n0+n1
 d=50

 dstar=3  # maximum number of dimensions allowable for our low-dimensional subspace
 Lhattarget=0.2  # maximum allowable loss for candidate low-dimensional subspace
 eq7w=1/2  # w from equation (7)
 rhostar=0.99 # rhostar is target posterior for morph
 eq6ws=rep(1,d)  # ws from equation (6) # multiplicative cost for movement in each dimension
 eq6ds=rep(1,d)  # ds from equation (6) # additive cost of doing *anything* in each dimension

gendat = function(myseed=12345)
## Roman suggested mixtures for each class ...
## so class 0 has two components and class 1 has five components.
## (1) the class 0's are rare, and
## (2) two of the dimensions (the first two, wlog) have signal.
#### (this probably requires n0 to be a multiple of 2 and n1 to be a multiple of 15)
 {
 # n0,n1,n,d are global variables
 require(MASS)
 set.seed(myseed)
 y=c(rep(0,n0),rep(1,n1))
  n01=n02=n0/2
  mu01=c(-5,-5,rep(0,d-2))
  mu02=c(+2,+2,rep(0,d-2))
  sigma01=sigma02=diag(c(2,2,rep(1,d-2)))
  x01 = mvrnorm(n01,mu01,sigma01)
  x02 = mvrnorm(n02,mu02,sigma02)
  n11=1*n1/15
  n12=2*n1/15
  n13=3*n1/15
  n14=4*n1/15
  n15=5*n1/15
  mu11=c(-1,-8,rep(0,d-2))
  mu12=c(-1,+4,rep(0,d-2))
  mu13=c(+0,+6,rep(0,d-2))
  mu14=c(+7,-1,rep(0,d-2))
  mu15=c(-3,-2,rep(0,d-2))
  sigma11=sigma12=sigma13=sigma14=sigma15=diag(c(2,2,rep(1,d-2)))
  x11 = mvrnorm(n11,mu11,sigma11)
  x12 = mvrnorm(n12,mu12,sigma12)
  x13 = mvrnorm(n13,mu13,sigma13)
  x14 = mvrnorm(n14,mu14,sigma14)
  x15 = mvrnorm(n15,mu15,sigma15)
  x = rbind(x01,x02,x11,x12,x13,x14,x15)
 return(list(x,y))
 }

 dat = gendat()
 x=dat[[1]]
 y=dat[[2]]
  plot(x[,1:2],col=y+1)  # scatterplot shows signal dimensions ...
 ## plot(x[,3:4],col=y+1)  # ... and (example) noise dimensions.

 #### make some plots ...
   ## mclust0 = Mclust(as.matrix(x[y==0,c(1,2)]) ,modelNames=c("EII","VII"))
   ## mclust1 = Mclust(as.matrix(x[y==1,c(1,2)]) ,modelNames=c("EII","VII"))

   ## plot(mclust0 , as.matrix(x[y==0,c(1,2)]) , what = "BIC")
   ## plot(mclust0 , as.matrix(x[y==0,c(1,2)]) , what = "classification")
   ## #plot(mclust0 , as.matrix(x[y==0,c(1,2)]) , what = "uncertainty")
   ## plot(mclust0 , as.matrix(x[y==0,c(1,2)]) , what = "density")
   ## points(-1.4,-2.05,cex=9,pch=".",col=2) # ????

   ## plot(mclust1 , as.matrix(x[y==1,c(1,2)]) , what = "BIC")
   ## plot(mclust1 , as.matrix(x[y==1,c(1,2)]) , what = "classification")
   ## #plot(mclust1 , as.matrix(x[y==1,c(1,2)]) , what = "uncertainty")
   ## plot(mclust1 , as.matrix(x[y==1,c(1,2)]) , what = "density")
   ## points(-1.91,-1.36,cex=9,pch=".",col=2) # ????
 ####

LhatResub = function(dat,w=eq7w)
 # w>pi0 up-weights loss impact for the class 0 (rare) cases;
 # w=1/2 yields a loss function which treats *percentage* misclassifications the same.
 {
 # n0,n1,n,d are global variables
 require(mclust)
 x=dat[[1]]
 y=dat[[2]]
  if(ncol(x)==1)
   {
   mclust0 = Mclust(as.matrix(x[y==0,]))
   mclust1 = Mclust(as.matrix(x[y==1,]))
   }
  if(ncol(x)>1)
   {
   mclust0 = Mclust(as.matrix(x[y==0,]),modelNames=c("EII","VII"))
   mclust1 = Mclust(as.matrix(x[y==1,]),modelNames=c("EII","VII"))
   }
 thisclass=NULL
 for(i in 1:nrow(x))
  {
  if(ncol(x)==1)
   {
   dens0 = dens(mclust0$modelName,x[i,],parameters=mclust0$parameters)
   dens1 = dens(mclust1$modelName,x[i,],parameters=mclust1$parameters)
   }
  if(ncol(x)>1)
   {
   dens0 = dens(mclust0$modelName,t(as.matrix(x[i,])),parameters=mclust0$parameters)
   dens1 = dens(mclust1$modelName,t(as.matrix(x[i,])),parameters=mclust1$parameters)
   }
  thisclass[i] = 1*((n1/n)*dens1 > (n0/n)*dens0) # i'm not worrying about loo priors.
  }
 Lhat = w*(1/n0)*sum(thisclass[y==0])+(1-w)*(1/n1)*sum(1-thisclass[y==1])
 return(Lhat)
 }

LhatLOO = function(dat,w=eq7w)
 # w>pi0 up-weights loss impact for the class 0 (rare) cases;
 # w=1/2 yields a loss function which treats *percentage* misclassifications the same.
 {
 # n0,n1,n,d are global variables
 require(mclust)
 x=dat[[1]]
 y=dat[[2]]
 thisclass=NULL
 for(i in 1:nrow(x))
  {
  xx = as.matrix(x[-i,])
  yy = y[-i]
  if(ncol(xx)==1)
   {
   mclust0 = Mclust(as.matrix(xx[yy==0,]))
   mclust1 = Mclust(as.matrix(xx[yy==1,]))
   dens0 = dens(mclust0$modelName,x[i,],parameters=mclust0$parameters)
   dens1 = dens(mclust1$modelName,x[i,],parameters=mclust1$parameters)
   }
  if(ncol(xx)>1)
   {
   mclust0 = Mclust(as.matrix(xx[yy==0,]),modelNames=c("EII","VII"))
   mclust1 = Mclust(as.matrix(xx[yy==1,]),modelNames=c("EII","VII"))
   dens0 = dens(mclust0$modelName,t(as.matrix(x[i,])),parameters=mclust0$parameters)
   dens1 = dens(mclust1$modelName,t(as.matrix(x[i,])),parameters=mclust1$parameters)
   }
  thisclass[i] = 1*((n1/n)*dens1 > (n0/n)*dens0) # i'm not worrying about loo priors.
  }
 Lhat = w*(1/n0)*sum(thisclass[y==0])+(1-w)*(1/n1)*sum(1-thisclass[y==1])
 return(Lhat)
 }

 ## Marginal Analysis: an approach a la [DGL p 566]:
 ## (we do this marginal analysis via *resubstitution* (for computational reasons).)
RL1d=NULL
#L1d=NULL
for(i in 1:d)
{
datp = list( as.matrix(dat[[1]][,i]) , dat[[2]] )
RL1d[i]=LhatResub(datp)
print(paste(i,round(RL1d[i],4)))
#L1d[i]=LhatLOO(datp)
#print(paste(i,round(L1d[i],4)))
}

dvec = order(RL1d)[1:dstar]
Lvec = sort(RL1d)[1:dstar]
 # > dvec
 # [1]  2 12 1
 # > Lvec
 # [1] 0.3850 0.4317 0.4533
 ## 'dvec' orders the canonical dimensions
 ## based on loss for mixture classifier.
 ## So (correctly, for this case) dims 1 & 2 are among the dstar=3 best;
 ## dim 12 is also under consideration.
 ## NB: None of the univariate marginals give low error.

 ## Joint Analysis:
 ## Now, consider the joint space given by dprime of these,
 ## for dprime=1,2,...,dstar.
 ## (The dprime=1 cases are done above, as 1-dimensional marginals.)
 ## (we do this joint analysis via *leave-one-out*.)

 # the full 50-dimensional space gives
LhatLOO(dat)
 # [1] 0.2667
 ## NB: full space does not give low error.

 # d'=1
for(i1 in 1:dstar)
{
datp = list( as.matrix(dat[[1]][,dvec[c(i1)]]) , dat[[2]] )
LLL=LhatLOO(datp)
print(paste(dvec[i1],round(LLL,4)))
}
 # [1] " 2 0.4133"
 # [1] "12 0.5067"
 # [1] " 1 0.4783"
 # (does not agree with original marginal search results above: this is loo vs resub.)

 # d'=2
for(i1 in 1:(dstar-1))
for(i2 in (i1+1):dstar)
{
datp = list( dat[[1]][,dvec[c(i1,i2)]] , dat[[2]] )
LLL=LhatLOO(datp)
print(paste(dvec[i1],dvec[i2],round(LLL,4)))
}
 # [1] " 2 12 0.3550"
 # [1] " 2  1 0.1667"   # <<<---!!! :-)
 # [1] "12  1 0.4700"
 # the pair 1 & 2, which should be the best, is the best.

 # d'=3
for(i1 in 1:(dstar-2))
for(i2 in (i1+1):(dstar-1))
for(i3 in (i2+1):dstar)
{
datp = list( dat[[1]][,dvec[c(i1,i2,i3)]] , dat[[2]] )
LLL=LhatLOO(datp)
print(paste(dvec[i1],dvec[i2],dvec[i3],round(LLL,4)))
}
 # [1] " 2 12 1 0.1667"   # <<<---!!! :-)

 # We need to have confidence that the low-dimensional subspace
 # admits a classifier with which we have confidence.
 # No 1-d subspace does.
 # The 2-d subspace given by dimensions (2,1) does.
 # So does the 3-d subspace given by dimensions (2,12,1).
 # So it is here that we will try to morph.

morph = function( myd , myx , rs=rhostar , dcost=rep(1,d) , dweight=rep(1,d) )
# rs is target posterior for morph
# dcost is additive cost of doing *anything* in each dimension
# dweight is multiplicative cost for movement in each dimension
 {
   myx = myx[myd]
   xx = dat[[1]]
   yy = dat[[2]]
   mclust0 = Mclust(as.matrix(xx[yy==0,myd]),modelNames=c("EII","VII"))
   mclust1 = Mclust(as.matrix(xx[yy==1,myd]),modelNames=c("EII","VII"))
   dens0 = dens(mclust0$modelName,t(as.matrix(myx)),parameters=mclust0$parameters)
   dens1 = dens(mclust1$modelName,t(as.matrix(myx)),parameters=mclust1$parameters)
   thispost = dens1/(dens0+dens1)
   minindexsofar=0
   mindistsofar=Inf
   for(i in 1:n1)
    {
    thisx = xx[n0+i,myd]
    dens0 = dens(mclust0$modelName,t(as.matrix(thisx)),parameters=mclust0$parameters)
    dens1 = dens(mclust1$modelName,t(as.matrix(thisx)),parameters=mclust1$parameters)
    post = dens1/(dens0+dens1)
    distsq = sum( dweight[myd]*(myx-thisx)^2 ) # using dist^2 ... just because
    if(post > rs) if(distsq<mindistsofar){ minindexsofar=i ; mindistsofar=distsq }
    }
   xp = xx[n0+minindexsofar,myd]
   mcost = sqrt(mindistsofar) + sum(dcost[myd])
   return(list(myd , myx , thispost,minindexsofar,mindistsofar, xp,mcost))
 }

 myx = c( c(-5,-5) , rep(0,d-2) )

for(i1 in 1:(dstar-1))
for(i2 in (i1+1):dstar)
 {
 datp = list( dat[[1]][,dvec[c(i1,i2)]] , dat[[2]] )
 LLL=LhatResub(datp)   # <<<--- resub!
 if(LLL <= Lhattarget)
  {
  print(paste(dvec[i1],dvec[i2],round(LLL,4)))
  xxx = morph( dvec[c(i1,i2)] , myx )
  print(xxx)
  }
 }
 # [1] "2 1 0.1383"
 # [[1]] [1] 2 1
 # [[2]] [1] -5 -5
 # [[3]] [1] 0.009747364
 # [[4]] [1] 3
 # [[5]] [1] 21.68562
 # [[6]] [1] -8.729351 -2.211171
 # [[7]] [1] 6.657

datp = list( dat[[1]][,dvec[c(1,2,3)]] , dat[[2]] )
LhatResub(datp)   # <<<--- resub!
# [1] 0.1166667
morph( dvec[c(1,2,3)] , myx )
 # [[1]] [1]  2 12  1
 # [[2]] [1] -5  0 -5
 # [[3]] [1] 0.004169471
 # [[4]] [1] 112
 # [[5]] [1] 15.97327
 # [[6]] [1] -1.7867857 -0.1814763 -2.6302755
 # [[7]] [1] 6.997

   # Given myx, we consider x' = arg min morphcost(x,x')
   # in spaces that meet the Lhattarget=0.2.
   # (That's just spaces (2,1) and (2,12,1).)
   # Of the class 1 observations that meet the posterior target,
   #  [[6]] [1] -8.729351 -2.211171
   # in dimensions
   #  [[1]] [1]  2 1
   # is least cost morph (6.657).
   # (Even though the 3-d morph is less distance (3.997 < 4.657)
   # it is more dimensions and this yields more total morph cost.)