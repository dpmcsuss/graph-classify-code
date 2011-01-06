mysimulator <- function(psi=NULL, init=NULL, LHS=NULL, RHS=NULL){       

  sigma <- function(y) {
    1
  }

  mu <- function(y) {
    cur.vol = psi[1]
    cur.beta = psi[2]
    cur.pi = psi[3]
    
    if( y*cur.vol > pi/2 ) {
      x = 0.999
    } else if (y*cur.vol < -pi/2) {
      x = 0.001
    } else {
      x = (sin(cur.vol*y)+1)/2 
    }
    
    L = (cur.beta/cur.vol)*(x-cur.pi)
    R = (-1/4)*cur.vol*(1-2*x)
    D = sqrt(x-x^2)
    
    (L+R)/D
  }
  
  dT = c(0,rep((RHS-LHS)/nstep,nstep))
  time.grid = LHS + cumsum(dT)
  Y = NULL;Y[1] = asin(2*init-1)/psi[1];
  dY = NULL;dY[1] = 0;

  for (i in 2:(length(time.grid))) {
    drift = mu(Y[i-1])*(1/nstep)
    diffs = rnorm(1)*sqrt(1/nstep)*sigma(Y[i-1])
    
    dY[i] = drift + diffs
    Y[i] = Y[i-1] + dY[i]   
  }

  X = (sin(psi[1]*Y)+1)/2

  retObj = list(LHS=LHS, RHS=RHS, start=init, end=X[[length(time.grid)]], ts.time=time.grid, ts.state=X)  
}

mysimulator.cond <- function(ts.obs.ij=NULL,cur.psi=NULL,maxTime=NULL,myunitRate=NULL) {

  lhs = 0
  rhs = ts.obs.ij[1,1]
  curinit.i = cur.psi[1,3]
  curinit.j = cur.psi[2,3]

  cur.init = numeric(nvertex)
  for(itr.v in 1:nvertex) {
    cur.init[itr.v] = cur.psi[itr.v,3]
  }

  times.i = NULL
  times.j = NULL
  states.i = NULL
  states.j = NULL
  
  for(myitr in 1:length(ts.obs.ij[,1])) {
      MORE = TRUE
      while(MORE) {
           
          cur.path.i = mysimulator(cur.psi[1,], init=curinit.i, LHS=lhs, RHS=rhs)
          cur.path.j = mysimulator(cur.psi[2,], init=curinit.j, LHS=lhs, RHS=rhs)

          
          temp.path.i = cbind(cur.path.i$ts.state,1-cur.path.i$ts.state,deparse.level=0)
          temp.path.j = cbind(cur.path.j$ts.state,1-cur.path.j$ts.state,deparse.level=0)
          
          temp.lambda.ij.first = as.vector(apply(cbind(temp.path.i[,1],temp.path.j[,1]),1,prod))
          temp.lambda.ij.second = as.vector(apply(cbind(temp.path.i[,2],temp.path.j[,2]),1,prod))
          temp.lambda.ij.both = cbind(temp.lambda.ij.first,temp.lambda.ij.second,deparse.level=0)
          temp.lambda.ij = as.vector(apply(temp.lambda.ij.both,1,sum))
          temp.end.i = c(cur.path.i$end,1-cur.path.i$end)
          temp.end.j = c(cur.path.j$end,1-cur.path.j$end)

          temp.nrow = nrow(temp.lambda.ij.both)
          temp.topic = ts.obs.ij[myitr,2]
          temp.actors = ts.TKVs[myitr,-c(1,2)]
          
          temp.lambda.ij.end = temp.lambda.ij.both[temp.nrow,temp.topic]

              
          ACCEPT = (runif(1) < temp.lambda.ij.end * exp(-myunitRate*mean(temp.lambda.ij)*(rhs-lhs)))

                  if(ACCEPT) {
#                    cat("******ACCEPT\n")
#                    cat("Topic ",temp.topic,"\n")
#                    cat("Both ",temp.lambda.ij.both[temp.nrow,],"\n")
            MORE = FALSE
          } else {
#                    cat("REJECT\n")
                    MORE = TRUE
                  }
      }

          
          if( myitr < length(ts.obs.ij[,1]) ) {
            lhs = rhs
            rhs = ts.obs.ij[myitr+1,1]
            curinit.i = cur.path.i$end
            curinit.j = cur.path.j$end
          } else {
            lhs = rhs
            rhs = maxTime
            curinit.i = cur.path.i$end
            curinit.j = cur.path.j$end
          }
          
          times.i = c(times.i,(cur.path.i$ts.time)[-1])
          times.j = c(times.j,(cur.path.j$ts.time)[-1])
          states.i = c(states.i,(cur.path.i$ts.state)[-1])
          states.j = c(states.j,(cur.path.j$ts.state)[-1])
  }

  MORE = TRUE
  while(MORE) {
    
    cur.path.i = mysimulator(cur.psi[1,], init=curinit.i, LHS=lhs, RHS=rhs)
    cur.path.j = mysimulator(cur.psi[2,], init=curinit.j, LHS=lhs, RHS=rhs)
    cur.path.lambda.ij = as.vector(apply(cbind(cur.path.i$ts.state, cur.path.j$ts.state),1,prod))

    temp.path.i = cbind(cur.path.i$ts.state,1-cur.path.i$ts.state,deparse.level=0)
    temp.path.j = cbind(cur.path.j$ts.state,1-cur.path.j$ts.state,deparse.level=0)
    
    temp.lambda.ij.first = as.vector(apply(cbind(temp.path.i[,1],temp.path.j[,1]),1,prod))
    temp.lambda.ij.second = as.vector(apply(cbind(temp.path.i[,2],temp.path.j[,2]),1,prod))
    temp.lambda.ij.both = cbind(temp.lambda.ij.first,temp.lambda.ij.second,deparse.level=0)    
    temp.lambda.ij = as.vector(apply(temp.lambda.ij.both,1,sum))
    
    ACCEPT = (runif(1) < exp(-myunitRate*mean(temp.lambda.ij)*(rhs-lhs)))
    
    if(ACCEPT) {
#      cat("******ACCEPT\n")
      MORE = FALSE
    } else {
#      cat("REJECT\n")
      MORE = TRUE
    }
  }

  cat("****************************FINISHED\n")
  times.i = c(times.i,cur.path.i$ts.time)
  times.j = c(times.j,cur.path.j$ts.time)
  states.i = c(states.i,cur.path.i$ts.state)
  states.j = c(states.j,cur.path.j$ts.state)

#  plot(times.i,states.i,type="l")
#  plot(times.j,states.j,type="l")
  
  retObj = list(ts.i=cbind(times.i,states.i), ts.j=cbind(times.j,states.j))
}

mydensity.vertex <-function(psi,grid.time,path) {
    first.index = 1
    last.index = length(path)    

    
#### HENG
    A = -(psi[2]/psi[1]^2+1/2)*(log(2) + log(path-path^2))
    B = psi[2]*(1-2*psi[3])/(2*psi[1]^2)*log(path/(1-path))
    C1 = A + B

    A= (psi[2]/psi[1] + psi[1]/2) * (2*path - 1)/(2*sqrt(path-path^2))
    B= (psi[2]*(1-2*psi[3]))/psi[1] * 1/(2*sqrt(path-path^2))
    C2 = (A + B)^2

    A = (psi[2]+psi[1]^2/2)/(4*(path-path^2))
    B = psi[2]*(1-2*psi[3])/psi[1]*(2*path-1)/(4*(path-path^2))
    C3 = A + B
######################


#
################## NAM
#    A = (psi[2]/psi[1]^2+1/2)*log(1/(2*sqrt(path-path^2)))
#    B = psi[2]*(-1*psi[3])/(psi[1]^2)*log(path/(1-path))
#    C1 = A + B
#
#    A= (psi[2]/psi[1]) * (path - psi[3])/sqrt(path-path^2)
#    B= (-1/2)*(psi[1]/2) * (1-2*path)/(2*sqrt(path-path^2))
#    C2 = (A + B)^2
#    
#    A = (psi[2]+psi[1]^2/2)/(4*(path-path^2))
#    B = psi[2]*(-2*psi[3])*(2*path-1)/(4*(path-path^2))
#    C3 = A + B
#################
#
#
    dretVal = C1[last.index] - C1[first.index]
    retVal = dretVal
#    
    dretVal = (-0.5)*sum(C2*c(0, diff(grid.time)))
    retVal = retVal + dretVal
#
    dretVal = (-0.5)*sum(C3*c(0, diff(grid.time)))
    retVal = retVal + dretVal


###############
#  

    path = (1/psi[1])*asin(2*path-1)

    A = (psi[2]/psi[1]^2+1/2)*log(1/(cos(psi[1]*path)))
    B = psi[2]/psi[1]^2*(1-2*psi[3])*log(1/cos(psi[1]*path)+tan(psi[1]*path))
    C1x = A + B

    A= (psi[2]/psi[1]+psi[1]/2) * tan(psi[1]*path)
    B= psi[2]/psi[1]*(1-2*psi[3])* 1/cos(psi[1]*path)
    C2x = (A + B)^2
    
    A = (psi[2]+psi[1]^2/2)/(cos(psi[1]*path))^2
    B = psi[2]/psi[1]*(1-2*psi[3])*(1/(cos(psi[1]*path)))*tan(psi[1]*path)
    C3x = A + B
#################

    retValx = C1x[last.index] - C1x[first.index]
    
    dretValx = (-0.5)*sum(C2x*c(0, diff(grid.time)))
    retValx = retValx + dretValx
    dretValx = (-0.5)*sum(C3x*c(0, diff(grid.time)))
    retValx = retValx + dretValx

#
#    cat(dretVal,"from orginal \n")
#    cat(retValx,"from L \n")
#    cat(dretVal,"from orginal \n")
#    cat(dretValx,"from L \n")
#    cat(dretVal,"from orginal \n")        
#    cat(dretValx,"from L \n\n")    



    as.numeric(retValx)
}
