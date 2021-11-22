library(igraph)
##Dataset from https://snap.stanford.edu/data/higgs-twitter.html 
##follower network "higgs-social_network.edgelist"
higgs.social_network <- read.table("~/networksampling/higgs-social_network.edgelist", 
                                   quote="\"", comment.char="")
##g is the raw data read and reserved in Rstudio work space
g<- graph.edgelist(as.matrix(higgs.social_network),directed=T)
A<-as_adjacency_matrix(g)## adjacency matrix
##Best sample size tested by experiments --15% of the whole agents. 
b<- 62072

##MHRW
mhrw_graph<- MHRW(g,b)
#Random Walk with Jump
rwwj_graph<- RWwJ_Sampling(g,b)
#Snowball Sampling
sbs_graph<-SBS(g,b)

##calculate weight
##for MHRW, targeting distribution is uniform distribution, there is no need to calculate 
##weights

##for random walk
##calculating weights, which are inverse inclusion probabilities
Z<-sum(1/(degree(rwwj_graph,mode = 'out')+0.3))
weight_rwwj<- (1/(degree(rwwj_graph)+0.3)/Z)

##for Snowball Sampling
####inclusion probability
N<-dim(A)[1]
denominator <- choose(N,3)
A2<-A%*%A
inclusion_prob<- c()

for (i in sbs_sample) {
  s1<- which(A[i,]==1)
  s2<-which(A2[i,]==1)
  ancestor<- union(s1,s2)
  numerator<-choose(N-length(ancestor),3)
  inclusion_prob<-append(1-numerator/denominator,inclusion_prob,after = 0)
}
weight_sbs<- 1/sbs.inclusion

###Estimation 

#1. MHRW

##in-degree distribution
M.in<-degree(mhrw_graph,mode = 'in')
M.dist.in <- degree_distribution(mhrw_graph, cumulative=F, mode="in")

##out-degree distribution
M.out<-degree(mhrw_graph,mode = 'out')
M.dist.out <- degree_distribution(mhrw_graph, cumulative=F, mode="out")

##Ratio
M.ratio<- M.in/M.out
M.ratio.dist<- table(M.ratio)/length(M.ratio)
M.N<-as.numeric(names(M.ratio.dist))
M.ratio.dist2<-as.data.frame(M.ratio.dist)
M.new<-M.ratio.dist2$Freq
for (i in 2:length(M.ratio.dist)) {
  M.new<-replace(M.new,i,M.new[i-1]+M.new[i])
}

##order
M.order<- length(V(mhrw_graph))


##mutual relationship
mutual <-which_mutual(mhrw_graph)
MUTUAL<-sum(mutual)
M.mutual<- (N/n)*MUTUAL

#2. RWwJ

R.number.in<-c()

for (i in 0: max(degree(rwwj_graph, mode = 'in'))) {
  vec<- which(degree(rwwj_graph, mode = 'in')==i)
  count<-sum(weight_rwwj[vec])*length(vec)
  R.number.in<- append(count,R.number.in,after = 0)
}

##number.in is the estimated number of nodes with in degree k
##R.frac.in is the estimated fraction
R.frac.in<-R.number.in/length(V(rwwj_graph))


##estimated out degree
R.number.out<- c()

for (j in 0:max(degree(rwwj_graph,mode = 'out'))) {
  vec<- which(degree(rwwj_graph, mode = 'out')==j)
  count.out<-sum(weight_rwwj[vec])*length(vec)
  R.number.out<- append(count.out,R.number.out,after = 0)
}
##number.out is the estimated number of nodes with out degree k
##R.frac.out is the estimated fraction
R.frac.out<- R.number.out/length(V(rwwj_graph))

##estimated ratio
R.ratio<-degree(rwwj_graph,mode = 'out')/degree(rwwj_graph,mode = 'in')
R.number.r<- c()

R.ratio<- R.ratio[which(is.finite(R.ratio))]
R.ratio<- na.omit(R.ratio)
sort(R.ratio)

for (k in R.ratio) {
  vec<-which(R.ratio==k)
  count.r<-sum(weight_rwwj[vec])*length(vec)
  R.number.r<- append(count.r,R.number.r,after = 0)
}
R.ratio.frac<-R.number.r/length(V(rwwj_graph))

##order
R.estimatedOrder<- sum(weight_rwwj)*length(V(rwwj_graph))

##mutual relationship
N_mutual<-which_mutual(rwwj_graph)
sum(N_mutual)

#3. SBS 
##estimated in degree 
S.number.in<-c()

for (i in 0:max(degree(sbs_graph,mode = 'in'))) {
  vec<- which(degree(sbs_graph, mode = 'in')==i)
  count<-sum(weight_sbs[vec])*length(vec)
  S.number.in<- append(count,S.number.in,after = 0)
}
##S.number.in is the estimated number of nodes with in degree k

S.frac.in<-S.number.in/length(V(sbs_graph))##fraction

S.number.out<- c()
for (j in 0:max(degree(sbs_graph,mode = 'out'))) {
  vec<- which(degree(sbs_graph, mode = 'out')==j)
  count.out<-sum(weight_sbs[vec])*length(vec)
  S.number.out<- append(count.out,S.number.out,after = 0)
}
##number.out is the estimated number of nodes with out degree k
S.frac.out<- S.number.out/length(V(sbs_graph))##fraction

##estimated ratio
S.ratio<-degree(sbs_graph,mode = 'out')/degree(sbs_graph,mode = 'in')
S.number.r<- c()

for (k in min(S.ratio):max(S.ratio)) {
  vec<-which(S.ratio==k)
  count.r<-sum(weight_sbs[vec])*length(vec)
  S.number.r<- append(count.r,S.number.r,after = 0)
}
S.ratio.frac<-S.number.r/length(V(sbs_graph))

##order
S.estimatedOrder<- sum(weight_sbs)*length(V(sbs_graph))

##mutual relationship
S.N_mutual<-which_mutual(sbs_graph)
sum(N_mutual)

###Comparison

###indegree
plot(x=0:max(deg.in), y=deg.dist.in, pch=4, cex=.5, col=1, 
     xlab="In-Degree", ylab="PDF",type = "p")
lines(x=0:max(M.in), y=M.dist.in, pch=2, cex=1.2, col=2,type = 'p')
lines(x=0:max(degree(rwwj_graph, mode = 'in')), y=R.frac.in,
      pch=19, cex=1.2, col=3,type = 'p')
lines(x=0:max(degree(sbs_graph, mode = 'in')),y=S.frac.in,pch=19, cex=1.2, col=4)
legend("topright",pch=c(15,15,15,15),legend=c("population","mhrw",'rwwj','sbs'),
       col=c(1,2,3,4),bty="n")

##outdegree
plot(x=0:max(deg.out), y=deg.dist.out, pch=19, cex=1.2, col=1, 
     xlab="Out-Degree", ylab="PDF",type = "l")
lines(x=0:max(M.out), y=M.dist.out, pch=19, cex=1.2, col=2)
plot(x=0:max(degree(rwwj_graph, mode = 'out')), y=R.frac.out, 
     pch=19, cex=1.2, col=3)
lines(x=0:max(degree(sbs_graph, mode = 'out')),y=S.frac.out,pch=19, cex=1.2, col=4)
legend("topright",pch=c(15,15,15,15),
       legend=c("population","mhrw",'rwwj','sbs'),col=c(1,2,3,4),bty="n")

#ratio
plot(x=as.numeric(names(population.ratio.dist)),
     y=population.new,ch=19, cex=1.2, col=1,
     xlab='follower vs. following ratio', ylab='distribution',type = "l")
lines(x=as.numeric(names(M.ratio.dist)),y=M.new, col=2,type = "l")
lines(x=R.ratio,y=R.ratio.frac, col=3,type = 'l')
lines(x=min(S.ratio):max(S.ratio) ,y=S.ratio.frac, col=4,type = 'l')
legend("topright",pch=c(15,15,15,15),
       legend=c("population","mhrw",'rwwj','sbs'),col=c(1,2,3,4),bty="n")


####EVALUATION
ks.test(deg.dist.in,M.dist.in)##in-degree evaluation
ks.test(deg.dist.out,M.dist.out)##out-degree evaluation

ks.test(deg.dist.in,R.frac.in)##in-degree evaluation
ks.test(deg.dist.out,R.frac.out)##out-degree evaluation

ks.test(deg.dist.in,S.frac.in)##in-degree evaluation
ks.test(deg.dist.out,S.frac.out)##out-degree evaluation

M.dist.in.new<-numeric(length(deg.dist.in))
M.dist.in.new[1]<-deg.dist.in[1]
for (i in 2:length(deg.dist.in)) {
  if(deg.in[i] %in% M.in){
    M.dist.in.new[i]<- deg.dist.in[i]
  }
  else M.dist.in.new[i]<-M.dist.in.new[i-1]
}

R.frac.in.new<-numeric(length(deg.dist.in))
R.frac.in.new[1]<-deg.dist.in[1]
for (i in 2:length(deg.dist.in)) {
  if(deg.in[i] %in% M.in){
    R.frac.in.new[i]<- deg.dist.in[i]
  }
  else R.frac.in.new[i]<-R.frac.in.new[i-1]
}

S.frac.in.new<-numeric(length(deg.dist.in))
S.frac.in.new[1]<-deg.dist.in[1]
for (i in 2:length(deg.dist.in)) {
  if(deg.in[i] %in% M.in){
    S.frac.in.new[i]<- deg.dist.in[i]
  }
  else S.frac.in.new[i]<-S.frac.in.new[i-1]
}

library(LaplacesDemon)
KLD(deg.dist.in,M.dist.in.new)$sum.KLD.py.px
KLD(deg.dist.in,R.frac.in.new)$sum.KLD.py.px
KLD(deg.dist.in,S.frac.in.new)$sum.KLD.py.px

M.dist.out.new<-numeric(length(deg.dist.out))
M.dist.out.new[1]<-deg.dist.out[1]
for (i in 2:length(deg.dist.out)) {
  if(deg.out[i] %in% M.out){
    M.dist.out.new[i]<- deg.dist.out[i]
  }
  else M.dist.out.new[i]<-M.dist.out.new[i-1]
}

R.frac.out.new<-numeric(length(deg.dist.out))
R.frac.out.new[1]<-deg.dist.out[1]
for (i in 2:length(deg.dist.out)) {
  if(deg.out[i] %in% R.out){
    R.frac.out.new[i]<- deg.dist.out[i]
  }
  else R.frac.out.new[i]<-R.frac.out.new[i-1]
}

S.frac.out.new<-numeric(length(deg.dist.out))
S.frac.out.new[1]<-deg.dist.out[1]
for (i in 2:length(deg.dist.out)) {
  if(deg.out[i] %in% M.out){
    S.frac.out.new[i]<- deg.dist.out[i]
  }
  else S.frac.out.new[i]<-S.frac.out.new[i-1]
}

library(LaplacesDemon)
KLD(deg.dist.out,M.dist.out.new)$sum.KLD.py.px
KLD(deg.dist.out,R.frac.out.new)$sum.KLD.py.px
KLD(deg.dist.out,S.frac.out.new)$sum.KLD.py.px


ks.test(R.ratio.frac, population.new)
ks.test(S.ratio.frac, population.new)
ks.test(M.ratio, population.new)


M.dist.in.new<-numeric(length(deg.dist.in))
M.dist.in.new[1]<-deg.dist.in[1]
for (i in 2:length(deg.dist.in)) {
  if(deg.in[i] %in% M.in){
    M.dist.in.new[i]<- deg.dist.in[i]
  }
  else M.dist.in.new[i]<-M.dist.in.new[i-1]
}

R.frac.in.new<-numeric(length(deg.dist.in))
R.frac.in.new[1]<-deg.dist.in[1]
for (i in 2:length(deg.dist.in)) {
  if(deg.in[i] %in% R.in){
    R.frac.in.new[i]<- deg.dist.in[i]
  }
  else R.frac.in.new[i]<-R.frac.in.new[i-1]
}

KLD(deg.dist.in,M.dist.in.new)$sum.KLD.py.px
KLD(deg.dist.in,R.frac.in.new)$sum.KLD.py.px
KLD(deg.dist.in,S.frac.in.new)$sum.KLD.py.px