
# X--hypergeometric with N units, M of them are good (N-M bad), select 
# n from N units, there are x good units in the sample.
# P(X=x)=dhyper(x,m,N-m,n)  for max(0,n+M-N)<=x <=min(n,M)
# Here, n and N are known, and want to estimate M based on x.
# i) one-sided CI for M, [L(X),N]
# ii) one-sided CI for M, [0,U(x)]
# iii) two-sided CI for M, [L(X),U(X)]
# iv) intervals for the proportion M/N can also be derived.
# phyper(x,a,b,n)=cdf of hypergeometric distribution with
# N=a+b units in the population, a good units, b bad units
# x good units in a sample of size n
 
# Let C1=[lcin1,ucin1] be the 1-alpha interval(=intersection of two one-sided
# 1-alpha/2 interval). We improve C1(=C^M_O in the paper), and let the resultant be
# Cw=[lciw,uciw](=C^M_I in the paper). Then C2 is a subset of Cw, and Cw is a subset of C1, where
# C2 is the 1-2alpha interval(=intersection of two one-sided
# 1-alpha interval).

tim=proc.time() 
N=400
alpha=0.05
n=200 # n must be an even number

# indicator fucntion of interval [a,b]
ind<- function(x,a,b){
(x>=a)*(x<=b)
}

lci<- function(x,n,alpha){   # the lower limit of 1-alpha CP interval
kk=1:length(x)
for (i in kk){
if (x[i]<0.5){kk[i]=0} else {aa=0:N
      bb=aa+1
      bb[2:(N+1)]=phyper(x[i]-1,aa[2:(N+1)]-1,N-aa[2:(N+1)]+1,n)
      dd=cbind(aa,bb)
      dd=dd[which(dd[,2]>=1-alpha),]
      if (length(dd)==2){kk[i]=dd[1]} else {kk[i]=max(dd[,1])}
      }}
kk
}

uci<- function(x,n,alpha){ # upper limit of 1-alpha CP interval
N-lci(n-x,n,alpha)
}

xx=0:n
lcin1=lci(xx,n,alpha/2)
ucin1=uci(xx,n,alpha/2)   #two-sided 1-alpha interval
lcin2=lci(xx,n,alpha)     
ucin2=uci(xx,n,alpha)   #two-sided 1-2alpha interval

lciw=lcin1  #  lcin1<=lciw<=lcin2
uciw=ucin1  #  ucin2<=uciw<=ucin1

# for an even number n 

xvalue=n/2+1  # start from the center
aa=lciw[xvalue]:floor(N/2)


ii=1
while (ii <length(aa)+0.5){ #ii
lciw[xvalue]=aa[ii]
uciw[xvalue]=N-aa[ii]

# the coverage probability function for the combined interval. 
cpci<- function(M){
kk=1:length(M)
for (i in kk){
xx<- 0:n
indp=xx
uu=0 
while (uu<n+0.5) {
indp[uu+1]=ind(M[i],lciw[uu+1],uciw[uu+1])*dhyper(uu,M[i],N-M[i],n)
uu=uu+1}
kk[i]=sum(indp)
}
kk
}
M=0:N
bb=min(cpci(M))
if (bb>=1-alpha){ii1=ii
ii=ii+1} else {ii=length(aa)+1}
} #ii
lciw[xvalue]=aa[ii1]
uciw[xvalue]=N-lciw[xvalue]


xvalue=n/2  # xvalue >=1
while (xvalue>0.5){  ##
al=lcin2[xvalue]-lciw[xvalue]+1
au=uciw[xvalue]-ucin2[xvalue]+1

if (al*au>1){ ### 

ff<-array(,dim=c(al*au,4))
for(i in 1:al)
{
ff[((i-1)*au+1):(i*au),1]=lciw[xvalue]+i-1
ff[((i-1)*au+1):(i*au),2]=(ucin2[xvalue]):(uciw[xvalue])
ff[((i-1)*au+1):(i*au),3]=ff[((i-1)*au+1):(i*au),2]-ff[((i-1)*au+1):(i*au),1]
}

for (ii in 1:dim(ff)[1]){ #ii
#for (ii in 1:1){ #ii

lciw[xvalue]=ff[ii,1]
uciw[xvalue]=ff[ii,2]
lciw[n+2-xvalue]=N-uciw[xvalue]
uciw[n+2-xvalue]=N-lciw[xvalue]

# the coverage probability function for the combined interval. 
cpci<- function(M){
kk=1:length(M)
for (i in kk){
xx<- 0:n
indp=xx
uu=0 
while (uu<n+0.5) {
indp[uu+1]=ind(M[i],lciw[uu+1],uciw[uu+1])*dhyper(uu,M[i],N-M[i],n)
uu=uu+1}
kk[i]=sum(indp)
}
kk
}
M=0:N
ff[ii,4]=min(cpci(M))

} #ii

ff=ff[which(ff[,4]>=1-alpha),]
if (length(ff)>4){
ff=ff[order(ff[,3]),]
lciw[xvalue]=ff[1,1]
uciw[xvalue]=ff[1,2]
} else {lciw[xvalue]=ff[1]
uciw[xvalue]=ff[2]}
lciw[n+2-xvalue]=N-uciw[xvalue]
uciw[n+2-xvalue]=N-lciw[xvalue]
} ###
xvalue=xvalue-1
}  ##
cbind(xx,lciw,uciw,lcin1,ucin1) # the improved and original intervals


(proc.time()-tim)/60   # computing time in minute 
