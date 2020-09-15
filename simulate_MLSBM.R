# Convert a membership matrix to a vector 

convertZ<-function(Z){
  n<-dim(Z)[1]
  return(sapply(1:n, function(i) {match(1, Z[i,])}))
}

# Convert a vector of labels to a membership matrix

convertClust<-function(clust){
  n<-length(clust)
  k<-length(unique(clust))
  Z<-matrix(0,nrow=n,ncol=k)
  for(i in 1:n){Z[i,clust[i]]<-1}
  return(Z)
}
# Generate a SBM
# n : number of nodes
# K : number of communities
# pi : KxK connectivity matrix
# convZ : community label vector


sampleSBM<-function(n, K, pi, convZ){
  simBern<-function(i,j){rbinom(1,1, p= pi[convZ[i],convZ[j]])}
  adjacency<-outer(1:n,1:n, Vectorize(simBern))
  adjacency[upper.tri(adjacency, diag=T)]<-0
  adjacency<-adjacency+t(adjacency)
  return(adjacency)
}


# Generate a MLSBM
# L : number of layers
# pi_l : tensor formed by the adjacency matrices of the different layers

sampleMLSBM<-function(n, K, pi_L, convZ){
  L<-dim(pi_L)[3]
  adj_tensor<-array(dim = c(n,n,L))
  for(l in 1:L){adj_tensor[,,l]<-sampleSBM(n,K,pi_L[,,l], convZ)}
  return(adj_tensor)
}

# Generate missing nodes
# rho : probability of selecting a node
# with_na : bolean indicating if the missing nodes should be filled with 0 or NA

## Generate missing nodes on a unilayer network
missSBM<-function(SBM_adj, rho, with_na=FALSE){
  obs_nodes<-rbinom(dim(SBM_adj)[1],1,rho)
  if(with_na==TRUE){obs_nodes<-replace(obs_nodes,obs_nodes==0,NA)}
  omega<-obs_nodes%o%obs_nodes
  SBM_adj_miss<-SBM_adj*omega
  liste<-list(SBM_adj_miss,obs_nodes)
  return(liste)
}

## Indicator function 
ind<-function(vec){
  sapply(vec,FUN=function(x)ifelse((x>=1),1,0))
}

## Generate missing nodes on a multilayer network
missMLSBM<-function(tensor_adj, rho, with_na=FALSE){
  L<-dim(tensor_adj)[3]
  n<-dim(tensor_adj)[1]
  tot_obs_nodes<-rep(0,n)
  obs_nodes_mat<-matrix(NA,nrow=n,ncol=L)
  for(l in 1:L){
    miss<-missSBM(tensor_adj[,,l],rho,with_na)
    tensor_adj[,,l]<-miss[[1]]
    obs_nodes_l<-miss[[2]]%>%replace_na(0)
    obs_nodes_mat[,l]<-obs_nodes_l
    tot_obs_nodes<-ind(tot_obs_nodes+obs_nodes_l)
  }
  return(list(tensor_adj,tot_obs_nodes,obs_nodes_mat))
}

# Delete the nodes that are not observed at least one time
delete_missing_nodes_tensor<-function(MLSBM_dataset,nodes){
  nodes<-as.logical(nodes)
  return(MLSBM_dataset[nodes,nodes,])
  
}

# Delete missing nodes in the labels
delete_missing_nodes_Z<-function(convZ,nodes){
  nodes<-as.logical(nodes)
  return(convZ[nodes])
  
}

# Generate connectivity matrices
# [a,b] is the range where the diagonal terms belong to
# r is such that the other terms are generating uniformly in [ra,rb]

conMat<-function(K,a,b,r){
  conMat<-matrix(runif(K^2,r*a,r*b),K,K)
  diag(conMat)<-runif(K,a,b)
  conMat<-forceSymmetric(conMat)
  return(conMat)
}

# Convert an array to a list
array_to_list<-function(MLSBM_dataset){
  n<-dim(MLSBM_dataset)[1]
  L<-dim(MLSBM_dataset)[3]
  list_adj<-vector("list",L)
  for(l in 1:L){
    list_adj[[l]]<-MLSBM_dataset[,,l]
  }
  return(list_adj)
}  