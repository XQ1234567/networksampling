#Random Walk with Jump

library(igraph)

RWwJ_Sampling<-function(G,b){
  A<-as_adjacency_matrix(G)
  s_0<- floor(runif(1,1,dim(A)[1]+1))
  v_i<- s_0
  k<-0
  alpha<- 300
  rwwj_sample<- c()
  edge_list<-c()
  while (k<b) {
    if (sum(A[v_i,])==0) {
      v_j<- floor(runif(1,1,dim(A)[1]+1))
      v_i<-v_j
      rwwj_sample<-append(v_j,rwwj_sample,after=0)
    }
    else{
      d<- rbinom(1,1,1-alpha/(sum(A[v_i,])+alpha))
      if (d==1) {
        d_v_i<- sum(A[v_i,])## out degree of node i
        new_vector<- which(A[v_i,]==1)
        v_j<-new_vector[floor(runif(1,1,d_v_i+1))]
        edge_list<-append(c(v_i,v_j),edge_list,after = 0)
        rwwj_sample<-append(v_j,rwwj_sample,after=0)
        v_i<-v_j
      }
      else{
        v_j<- floor(runif(1,1,dim(A)[1]+1))
        v_i<-v_j
        rwwj_sample<-append(v_j,rwwj_sample,after=0)
      }
    }
    k<-k+1
  }
  edge_list<-matrix(edge_list, ncol=2,byrow=TRUE)
  rwwj_graph<- graph.edgelist(edge_list)
  rwwj_graph<-simplify(rwwj_graph,remove.multiple = TRUE)
  differences<- setdiff(unique(rwwj_sample),V(rwwj_graph))
  rwwj_graph<- rwwj_graph%>%
    add_vertices(differences)
  return(rwwj_graph)
}

#rwwj_sample<-RWwJ_sampling(g,b)##sampling process finished
#rwwj_graph<- induced_subgraph(g,rwwj_sample)##rwwj_subgraph
#is.igraph(rwwj_graph)

