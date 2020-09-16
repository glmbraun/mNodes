###################################
# The MIT reaility mining dataset #
###################################

data(reality_mining_1392)

# Create a multilayer network by merging some layer (summation) and applying an indicator function
sum_layer<-function(week){
  list_mat<-array(0, dim=c(96,96))
  for(i in 1:42){
    list_mat<-list_mat+reality_mining_1392[,,24+42*(week-1)+i]
  }
  list_mat<-ifelse(list_mat>0, 1, 0)
  return(list_mat)
}


res<-lapply(1:32,sum_layer)

tensAdj<-array(NA,dim=c(96,96,32))
for(l in 1:32){
  tensAdj[,,l]<-res[[l]]
}

obs_nodes_mat<-matrix(1,nrow=96,ncol=32)

clust1<-clust_on_mat_eig(layer_sum,2)
clust2<-lmfo_miss(tensAdj,obs_nodes_mat,2)
clust3<-kpod(final_aggr_miss(tensAdj,2,obs_nodes_mat),2)$cluster
clust3<-kmeans(final_aggr_miss(tensAdj,2,obs_nodes_mat),2)$cluster

#kpod clust is very sensitive to initialization ? Error in final_aggr_miss : should use svds()
NMI(clust1,clust2)
NMI(clust2,clust3)

# What happens if we delete some nodes ?
# We use clust1 as the thrue underliying partition for convZ

simulMiss<-function(nb, rho, convZ, K=2){
  
  score<-matrix(NA,nrow=4,ncol=nb)
  rownames(score)<-c("sum_of_adj","k_pod","OLMF","imput")
  
  for(i in 1:nb){
    
    # Generate missing nodes
    miss<-missMLSBM(tensAdj,rho,with_na = TRUE)
    tensAdj_m<-miss[[1]]
    tot_obs_nodes<-miss[[2]]
    obs_nodes_mat<-miss[[3]]
    
    #  #Delete nodes that not have been observed at least one time
    tensorAdj_del<-delete_missing_nodes_tensor(tensAdj_m,tot_obs_nodes)
    convZ_del<-delete_missing_nodes_Z(convZ,tot_obs_nodes)
    obs_nodes_mat_del<-obs_nodes_mat[tot_obs_nodes%>%as.logical(),]
    tensorAdj_del0<-tensorAdj_del%>%replace_na(0)
    
    
    # Apply different clustering methods and compute the score
    ## Sum of adjacency matrices
    
    sum_of_adjMat<-sum_of_adj(tensorAdj_del0)
    clust<-clust_on_mat_eig(sum_of_adjMat,K)
    score[1,i]<-NMI(clust,convZ_del)
    
    ## kpod_clust
    # tot_obs_nodes<-tot_obs_nodes%>%as.logical()
    # final_ag<-final_aggr_miss(tensorAdj_del,K,obs_nodes_mat_del)
    # clust<-kpod(final_ag,K)$cluster
    # score[2,i]<-NMI(clust,convZ_del)
    
    ## OLMF modified for missing nodes
    n_obs<-dim(tensorAdj_del0)[1]
    #x<-array_to_list(tensorAdj_del0)
    #score[4,i]<-NMI(lmfo(x,n_obs,K),convZ_del)
    score[3,i]<-NMI(lmfo_miss(tensorAdj_del0,obs_nodes_mat_del,K),convZ_del)
    
    ##Imputation method
    tensorAdj_in<-tensorAdj_del0
    for(rep in 1:20){
      Z_in<-convertClust(clust)
      tensorAdj_out<-impute_nodes(Z_in,tensorAdj_in,obs_nodes_mat_del,K)
      sum_of_adjn<-sum_of_adj(tensorAdj_out)
      clust<-clust_on_mat_eig(sum_of_adjn,K)
      tensorAdj_in<-tensorAdj_out
    }
    score[4,i]<-NMI(clust,convZ_del)
    
  }
  return(t(rowMeans(score)))
}


score<-array(NA,dim=c(10,4))
score[1,]<-simulMiss(1,0.1,clust1)

"for(r in seq(0.1,1,by=0.1)){
  score[r*10,]<-simulMiss(50,r,clust1)
}

sapply(seq(0.1,1,by=0.1),function(x){return(simulMiss(2,x,clust1))})
