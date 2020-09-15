# Load packages
library(Matrix)
library(lbfgs)
library(ggplot2)
library(dplyr)
library(RSpectra)
library(aricode)
library(kpodclustr)
library(janitor)
library(tidyr)
library(purrr)
library(grid)
library(microbenchmark)
library(foreach)
library(doParallel)
library(igraph)
library(gridExtra)
library(lemon)
library(latex2exp)
library(GreedySBTM)

#############################
# Results on simulated data #
#############################

repSim_miss<-function(nb,K,nb_nodes,L,rho,a,b,r){
  score<-matrix(NA,nrow=6,ncol=nb)
  rownames(score)<-c("sum_of_adj","sum_of_sq","k_pod","OLMF","imput","tot_obs_nodes")
  pi_L<-array(NA,c(K,K,L))
  
  for(i in 1:nb){
    # Generate a MLSBM (with balanced communities)
    alpha<-rep(1/K,K)
    Z = t(rmultinom(nb_nodes,1, alpha))
    convZ<-convertZ(Z)
    for(j in 1:L){
      pi_L[,,j]<-array(conMat(K,a,b,r))
    }
    tensorAdj<-sampleMLSBM(nb_nodes,K,pi_L,convZ)
    
    # Generate missing nodes
    miss<-missMLSBM(tensorAdj,rho,with_na = TRUE)
    tensorAdj<-miss[[1]]
    tot_obs_nodes<-miss[[2]]
    obs_nodes_mat<-miss[[3]]
    score[6,i]<-sum(tot_obs_nodes)
    
    #Delete nodes that not have been observed at least one time
    tensorAdj_del<-delete_missing_nodes_tensor(tensorAdj,tot_obs_nodes)
    convZ_del<-delete_missing_nodes_Z(convZ,tot_obs_nodes)
    obs_nodes_mat_del<-obs_nodes_mat[tot_obs_nodes%>%as.logical(),]
    tensorAdj_del0<-tensorAdj_del%>%replace_na(0)
    
    # Apply different clustering methods and compute the score
    ## Sum of adjacency matrices
    
    sum_of_adjMat<-sum_of_adj(tensorAdj_del0)
    clust<-clust_on_mat_eig(sum_of_adjMat,K)
    score[1,i]<-NMI(clust,convZ_del)
    
    # Sum of squared adjacency matrices
    # sum_of_square<-sum_of_square_adj(tensorAdj_del0)
    # clust<-clust_on_mat_eig(sum_of_square,K)
    # score[2,i]<-NMI(clust,convZ_del)
    
    # kpod_clust
    tot_obs_nodes<-tot_obs_nodes%>%as.logical()
    final_ag<-final_aggr_miss(tensorAdj_del,K,obs_nodes_mat_del)
    clust<-kpod(final_ag,K)$cluster
    score[3,i]<-NMI(clust,convZ_del)
    
    # OLMF modified for missing nodes
    n_obs<-dim(tensorAdj_del0)[1]
    score[4,i]<-NMI(lmfo_miss(tensorAdj_del0,obs_nodes_mat_del,K),convZ_del)
    
    #Imputation method
    
    tensorAdj_in<-tensorAdj_del0
    for(rep in 1:20){
      Z_in<-convertClust(clust)
      tensorAdj_out<-impute_nodes(Z_in,tensorAdj_in,obs_nodes_mat_del,K)
      sum_of_adjn<-sum_of_adj(tensorAdj_out)
      clust<-clust_on_mat_eig(sum_of_adjn,K)
      tensorAdj_in<-tensorAdj_out
    }
    score[5,i]<-NMI(clust,convZ_del)
    
  }
  return(t(rowMeans(score)))
}

# test<-repSim_miss(1,3,600,2,0.8,0.18,0.19,0.5)

#############################
# Variying number of layers 
#############################

seq_l<-seq(3,6,by=3)

cl <- makeForkCluster(10)
registerDoParallel(cl)
result <- foreach(rho=1:10,.combine = rbind) %dopar% {
  lapply(seq_l,function(l){repSim_miss(1,3,1000,l,rho/10,0.18,0.19,0.7)%>%as.data.frame()%>%mutate(rho=rho/10,layer=l)})
}

stopCluster(cl)
registerDoSEQ()

data<-array(NA,dim=c(0,7))
for(i in 1:length(result)){
  data<-rbind(data,result[[i]])
}

###############################
# Variying number of nodes
###############################
# Experiment 2' : L=3, K= 3 fixed and n vary. The experiment is repeated for different values of rho.
# a= 0.18 b=0.19 r=0.7 so there is no overlap 

seq_n<-seq(600,2600,by=200)

cl <- makeForkCluster(10)
registerDoParallel(cl)
result2 <- foreach(rho=1:10,.combine = rbind) %dopar% {
  lapply(seq_n,function(n){repSim_miss(20,3,n,3,rho/10,0.18,0.19,0.7)%>%as.data.frame()%>%mutate(rho=rho/10,nodes=n)})
}
stopCluster(cl)
registerDoSEQ()

data2<-array(NA,dim=c(0,8))
for(i in 1:length(result2)){
  data2<-rbind(data2,result2[[i]])
}

#################################
# Plots
#################################

for(l in seq_l){
  assign(paste("d", l, sep = ""),data%>%filter(layer==l))
  
  pv<-ggplot(get(paste("d", l, sep = "")), aes(rho,label=tot_obs_nodes)) + 
    geom_line(aes(y = sum_of_adj, colour = "sum_of_adj",linetype="sum_of_adj")) + 
    geom_line(aes(y = k_pod, colour = "k_pod_clust",linetype = "k_pod_clust"))+
    geom_line(aes(y = OLMF, colour = "OLMF",linetype="OLMF"))+
    geom_line(aes(y = imput, colour = "Impute",linetype="Impute"))+
    scale_color_manual("",values=c("deepskyblue2","black","orange","deepskyblue3"))+
    scale_linetype_manual("",values=c("dotted","solid","longdash","solid"))+
    labs(y="NMI",x=TeX("$\\rho$"),title=paste("L=",l,", K=3, n=1000",sep=""))+
    
    theme(plot.margin = unit(c(1,1,1,1), "lines"),plot.title = element_text(hjust = 0.5, vjust=5)) +
    coord_cartesian(clip = "off")
  
  for(i in 1:10){
    pv<-pv + annotation_custom(
      grob = textGrob(label = round(get(paste("d", l, sep = ""))$tot_obs_nodes[i]), hjust = 0, gp = gpar(cex = 0.6)),
      xmin = get(paste("d", l, sep = ""))$rho[i],      # Vertical position of the textGrob
      xmax = get(paste("d", l, sep = ""))$rho[i],
      ymin = 1.05,         # Note: The grobs are positioned outside the plot area
      ymax = 1.05)
  }
  
  assign(paste("p",l,sep=""),pv)
  
}

grid_arrange_shared_legend(p3,p6,p9,p12,nrow=2,ncol=2)

## Varying number of layer for a given rho
seq_rho<-seq(0.1,1,by=0.1)
for(r in seq_rho){
  assign(paste("d", r, sep = ""),data%>%filter(rho==r))
  
  pv<-ggplot(get(paste("d", r, sep = "")), aes(layer,label=tot_obs_nodes)) + 
    geom_line(aes(y = sum_of_adj, colour = "sum_of_adj",linetype="sum_of_adj")) + 
    geom_line(aes(y = k_pod, colour = "k_pod_clust",linetype="k_pod_clust"))+
    geom_line(aes(y = OLMF, colour = "OLMF",linetype="OLMF"))+
    geom_line(aes(y = imput, colour = "Impute",linetype="Impute"))+
    scale_color_manual("",values=c("deepskyblue2","black","orange","deepskyblue3"))+
    scale_linetype_manual("",values=c("dotted","solid","longdash","solid"))+
    labs(y="NMI",x="L",title=TeX(sprintf("$\\rho = %1.1f$ K=3, n=1000",r)))+
    theme(plot.margin = unit(c(1,1,1,1), "lines"),legend.title=element_blank(),plot.title = element_text(hjust = 0.5, vjust=5)) +
    coord_cartesian(ylim=c(0,1),clip = "off")
  
  for(i in 1:4){
    pv<-pv + annotation_custom(
      grob = textGrob(label = round(get(paste("d", r, sep = ""))$tot_obs_nodes[i]), hjust = 0, gp = gpar(cex = 0.6)),
      xmin = get(paste("d", r, sep = ""))$layer[i],      # Vertical position of the textGrob
      xmax = get(paste("d", r, sep = ""))$layer[i],
      ymin = 1.05,         # Note: The grobs are positioned outside the plot area
      ymax = 1.05)
  }
  
  assign(paste("p",r,sep=""),pv)
  
}

grid_arrange_shared_legend(p0.4,p0.5,p0.6,p0.8,nrow=2,ncol=2)

# Variying number of nodes

for(r in seq_rho){
  assign(paste("d2", r, sep = ""),data2%>%filter(rho==r))
  
  pv<-ggplot(get(paste("d2", r, sep = "")), aes(nodes,label=tot_obs_nodes)) + 
    geom_line(aes(y = sum_of_adj, colour = "sum_of_adj",linetype="sum_of_adj")) + 
    geom_line(aes(y = k_pod, colour = "k_pod_clust",linetype = "k_pod_clust"))+
    geom_line(aes(y = OLMF, colour = "OLMF",linetype="OLMF"))+
    geom_line(aes(y = imput, colour = "Impute",linetype="Impute"))+
    scale_color_manual("",values=c("deepskyblue2","black","orange","deepskyblue3"))+
    scale_linetype_manual("",values=c("dotted","solid","longdash","solid"))+
    labs(y="NMI",x="n",title=TeX(sprintf("$\\rho = %1.1f$ K=3, n=1000",r)))+
    theme(plot.margin = unit(c(1,1,1,1), "lines"),legend.title=element_blank(),plot.title = element_text(hjust = 0.5, vjust=5)) +
    coord_cartesian(ylim=c(0,1),clip = "off")
  
  for(i in 1:10){
    pv<-pv + annotation_custom(
      grob = textGrob(label = round(get(paste("d2", r, sep = ""))$tot_obs_nodes[i]), hjust = 0, gp = gpar(cex = 0.6)),
      xmin = get(paste("d2", r, sep = ""))$nodes[i],      # Vertical position of the textGrob
      xmax =get(paste("d2", r, sep = ""))$nodes[i],
      ymin = 1.05,         # Note: The grobs are positioned outside the plot area
      ymax = 1.05)
  }
  
  assign(paste("pp",r,sep=""),pv)
  
}

grid_arrange_shared_legend(pp0.4,pp0.5,pp0.6,pp0.8,nrow=2,ncol=2)
