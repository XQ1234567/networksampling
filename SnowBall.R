##sampling process 3 Snowball Sampling

SBS<- function(G,b){
  sbs_sample <- c()
  seed<-3
  set.seed(100)## assume we have 100 seeds.
  starter <- sample(1:length(V(g)), seed)
  current <- c()
  current[1:seed] <- V(g)$name[starter]
  count <- seed
  sbs_sample[1:seed] <- current[1:seed]
  stage<- 0 
  while (stage < 3) {
    nnode = length(current)
    for (i in 1:nnode) {
      ngh <- neighbors(g, current[i],mode = "out")
      sbs_sample <- c(sbs_sample, V(g)$name[ngh])
      sbs_sample <- unique(sbs_sample)
    }
    tmp_sample <- sbs_sample[(count + 1):length(sbs_sample)]
    current <- tmp_sample
    count <- length(sbs_sample)
    stage<-stage+1
  }
  sbs_graph<- induced_subgraph(g,sbs_sample)
  return(sbs_graph)
}
