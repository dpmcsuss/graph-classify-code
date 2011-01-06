# re: Priebe's Conjecture #1
# this is assuming:
#  p and q are known, with q>p;
#  |\hatmcE|=k;
#  |\hatmcE \cap \mcE|=tt;
#  log((1-q)/(1-p))/log(p*(1-q)/q*(1-p)) is irrational.
# if my math is correct ...
L = function(tt,p=0.1,q=0.5,pi1=1/2)
 {
 x=log((1-q)/(1-p))/log(p*(1-q)/q*(1-p))
 L0 = (1-pbinom(x*tt,tt,p)) # P[g(G)=1|Y=0]
 L1 = pbinom(x*tt,tt,q)     # P[g(G)=0|Y=1]
 print(x*tt)
 print(L0)
 print(L1)
 print(pi1*L1 + (1-pi1)*L0) # P[g(G) \neq Y]
 }
# ... then L is *not* monotonic decreasing in tt!
# (thereby ***disproving*** Priebe's Conjecture #1.)
# to wit:
# > L(0:5)
# [1] 0.000 0.244 0.488 0.732 0.9764 1.22051
# [1] 0.001 0.100 0.190 0.271 0.3439 0.08146
# [1] 1.000 0.500 0.250 0.125 0.0625 0.18750
# [1] 0.500 0.300 0.220 0.198 0.2032 0.13448
# notice that L(tt=4)>L(tt=3)!?
# i'm not sure i understand ...
# i'm not sure i believe it ...
# i am sure i am soliciting help!
# ---cep
# 02/19/2010