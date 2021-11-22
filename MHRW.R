##Metropolis Hastings algorithm with Jumping

MHRW<-function(G,b){
  s_0<- floor(runif(1,1,dim(A)[1]))
  v_i<- s_0
  k<-0
  b<- 62072
  mhrw_sample<- c()
  while (k<b) {
    d_v_i<- sum(A[v_i,])## out degree of node i
    
    new_vector<- which(A[v_i,]==1)## which nodes are connected to v_i
    v_j<-new_vector[floor(runif(1,1,d_v_i))]
    d_v_j<- sum(A[v_j,])## out degree of node j
    p<- runif(1)
    accept_rate<- min(d_v_i/d_v_j,1)
    if (p<accept_rate) {
      v_i<-v_j
      mhrw_sample<-append(v_j,mhrw_sample,after=0)
    }
    else mhrw_sample<- append(v_i,mhrw_sample,after=0)
    k<-k+1
    
  }
  ##construct the subgraph
  mhrw_graph<- induced_subgraph(g,mhrw_sample)
  return(mhrw_graph)
}

