

require(igraph)

load("celegansGraph.Rd")
gg <- graph.adjacency(Ag,mode="undirected")
gg <- set.graph.attribute(gg,"animal","C. Elegans")
gg <- set.graph.attribute(gg,"type","Gap Junction")

gg <- set.vertex.attribute(gg,"category",value=c("M","I","S")[vcols])
gg <- set.vertex.attribute(gg,"type",value=types)
gg <- set.vertex.attribute(gg,"name",value=neurons)

ggu <- graph.adjacency(Ag>0,mode="undirected")
ggu <- set.graph.attribute(ggu,"animal","C. Elegans")
ggu <- set.graph.attribute(ggu,"type","Gap Junction")

ggu <- set.vertex.attribute(ggu,"category",value=c("M","I","S")[vcols])
ggu <- set.vertex.attribute(ggu,"type",value=types)
ggu <- set.vertex.attribute(ggu,"name",value=neurons)

gc <- graph.adjacency(Ac>0,mode="directed")
gc <- set.graph.attribute(gc,"animal","C. Elegans")
gc <- set.graph.attribute(gc,"type","Chemical Synapse")

gc <- set.vertex.attribute(gc,"category",value=c("M","I","S")[vcols])
gc <- set.vertex.attribute(gc,"type",value=types)
gc <- set.vertex.attribute(gc,"name",value=neurons)

gc <- set.edge.attribute(gc,"weight",value=Ac[get.edgelist(gc,names=F)+1])

spg <- shortest.paths(gg)
spg[is.infinite(spg)] <- 2+max(spg[!is.infinite(spg)])
spc <- shortest.paths(gc)
spc[is.infinite(spc)] <- 2+max(spc[!is.infinite(spc)])
spg <- spg/max(spg)
spc <- spc/max(spc)
M <- matrix(0,nrow=2*nrow(Ag),ncol=2*ncol(Ag))
M[1:279,1:279] <- spg
M[-(1:279),-(1:279)] <- spc
A <- (spg+spc)/2
M[(1:279),-(1:279)] <- A
M[-(1:279),(1:279)] <- t(A)
x <- cmdscale(M)
require(shapes)
z <- procOPA(x[1:279,],x[-(1:279),])
plot(z$Ahat,pch=20,xlim=range(c(z$Ahat,z$Bhat)),ylim=range(c(z$Ahat,z$Bhat)))
points(z$Bhat,pch=20,col=2)
segments(z$Ahat[,1],z$Ahat[,2],z$Bhat[,1],z$Bhat[,2],col=gray(0.7))

par(mfrow=c(2,2))
plot(z$Ahat,pch=20,cex=.5,xlab="",ylab="",main="Gap vs Chem",
     xlim=range(c(z$Ahat,z$Bhat)),ylim=range(c(z$Ahat,z$Bhat)))
points(z$Bhat,pch=20,col=2,cex=.5)
segments(z$Ahat[,1],z$Ahat[,2],z$Bhat[,1],z$Bhat[,2],col=gray(0.8))
points(z$Bhat,pch=20,col=2,cex=.5)
points(z$Ahat,pch=20,col=1,cex=.5)
legend(-.5,.6,legend=c("Gap","Chem"),col=1:2,lty=1)

plot(z$Ahat,pch=20,cex=.5,type="n",main="Gap vs Chem",
     xlim=range(c(z$Ahat,z$Bhat)),ylim=range(c(z$Ahat,z$Bhat)))
segments(z$Ahat[,1],z$Ahat[,2],z$Bhat[,1],z$Bhat[,2],col=gray(0.8))

plot(z$Ahat,pch=20,cex=.5,xlab="",ylab="",
     main=get.graph.attribute(gg,"type"),
     xlim=range(c(z$Ahat,z$Bhat)),ylim=range(c(z$Ahat,z$Bhat)))
a <- get.edgelist(gg,names=F)+1
segments(z$Ahat[a[,1],1],z$Ahat[a[,1],2],
         z$Ahat[a[,2],1],z$Ahat[a[,2],2],
			col=gray(.8))
points(z$Ahat,pch=20,col=1,cex=.5)

plot(z$Bhat,pch=20,cex=.5,xlab="",ylab="",
     main=get.graph.attribute(gc,"type"),
     xlim=range(c(z$Ahat,z$Bhat)),ylim=range(c(z$Ahat,z$Bhat)))
a <- get.edgelist(gc,names=F)+1
arrows(z$Bhat[a[,1],1],z$Bhat[a[,1],2],
         z$Bhat[a[,2],1],z$Bhat[a[,2],2],
			col=gray(.8),length=.05)
points(z$Bhat,pch=20,col=1,cex=.5)

#plot(gg,layout=z$Ahat,vertex.color=vcols,
#     vertex.label=get.vertex.attribute(gg,"type"))
#plot(gc,layout=z$Bhat,vertex.color=vcols,
#     vertex.label=get.vertex.attribute(gc,"type"))

