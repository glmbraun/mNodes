# Sum of adjacency matrices
sum_of_adj<-function(sampleMLSBM){
  sum_of_adj<-apply(sampleMLSBM,c(1,2),sum)
  return(sum_of_adj)
}

# Mean of observed entries
mean_of_obs_ent<-function(sampleMLSBM_na){
  n<-dim(sampleMLSBM_na)[1]
  fun_mean<-function(x,y){
    res<-sampleMLSBM_na[x,y,]%>%mean(,na.rm=TRUE)
    return(res)
  }
  sapply(1:n,function(i) sapply(1:n, function(j)fun_mean(i,j)))
}

sampleMLSBM_na<-missMLSBM(MLSBM_dataset,0.8,with_na=TRUE)[[1]]

# Clustering on a matrix
clust_on_mat_eig<-function(A,k){
  return(kmeans(svds(A,k)$u,k)$cluster)
}

# Sum of the square of adjacency matrices

sum_of_square_adj<-function(sampleMLSBM){
  L<-dim(sampleMLSBM)[3]
  n<-dim(sampleMLSBM)[1]
  sum_of_square<-matrix(0,nrow=n,ncol=n)
  for(l in 1:L){
    adj_mat<-as.matrix(sampleMLSBM[,,l])
    square<-adj_mat%*%adj_mat
    diag<-diag(square)
    square<-square-diag(diag)
    sum_of_square<-sum_of_square+square
  }
  return(sum_of_square)
}

# Aggregate spectral kernel

aggregate_kernel<-function(sampleMLSBM,k){
  aggregate_k<-matrix(0,n,n)
  for(l in 1:L){
    proj<-eigs_sym(sampleMLSBM[,,l],k)$vectors
    aggregate_k<-aggregate_k+proj%*%t(proj)
  }
  return(aggregate_k)
}



# k-means on the flat matrix (nxKL) obtained by stacking the eigenvectors of each layers

final_aggr_miss<-function(sampleMLSBM,k,obs_nodes_mat){
  n<-dim(sampleMLSBM)[1]
  L<-dim(sampleMLSBM)[3]
  obs_nodes_l<-obs_nodes_mat[,1]%>%as.logical()
  eig<-matrix(0,nrow=n,ncol=k)
  flat_eig<-matrix(NA,nrow=n,ncol=k)
  flat_eig[obs_nodes_l,]<-svd(sampleMLSBM[obs_nodes_l,obs_nodes_l,1],k)$u
  for(l in 2:L){
    obs_nodes_l<-obs_nodes_mat[,l]%>%as.logical()
    eig2<-matrix(NA,nrow=n,ncol=k)
    c<-sampleMLSBM[obs_nodes_l,obs_nodes_l,l]
    eig2[obs_nodes_l,]<-svd(sampleMLSBM[obs_nodes_l,obs_nodes_l,l],k)$u
    flat_eig<-cbind(flat_eig,eig2)
  }
  return(flat_eig)
}


miss<-missMLSBM(tensAdj,rho,with_na = TRUE)
tensAdj_m<-miss[[1]]
tot_obs_nodes<-miss[[2]]
obs_nodes_mat<-miss[[3]]

final_aggr_miss(tensorAdj_del,2,obs_nodes_mat_del )


# Missing nodes imputation

impute_nodes<-function(Z_in, tensorAdj,obs_nodes_mat,k){
  n<-dim(tensorAdj)[1]
  L<-dim(tensorAdj)[3]
  tensorAdj_imp<-array(dim = c(n,n,L))
  Lambda<-array(dim = c(k,k,L))
  J<-matrix(1,n,n)
  tensorAdj<-tensorAdj%>%replace_na(0)
  for(l in 1:L){
    #Pr<-Z_in%*%solve(t(Z_in)%*%Z_in)%*%t(Z_in)
    Lambda[,,l]<-(t(Z_in)%*%as.matrix(tensorAdj[,,l])%*%Z_in)/(t(Z_in)%*%J%*%Z_in)
    obs_nodes_l<-obs_nodes_mat[,l]%>%as.logical()
    tensorAdj_imp[,,l]<-Z_in%*%as.matrix(Lambda[,,l])%*%t(Z_in)
    tensorAdj_imp[obs_nodes_l,obs_nodes_l,l]<-tensorAdj[obs_nodes_l,obs_nodes_l,l]
  }
  return(tensorAdj_imp)
}

# sample<-missMLSBM(MLSBM_dataset, 0.7, TRUE)
# sample_tensor<-sample[[1]]
# sum_of<-sum_of_adj(sample_tensor%>%replace_na(0))
# clust<-clust_on_mat_eig(sum_of,k)
# Z_in<-convertClust(clust)
# obs_nodes_mat<-sample[[3]]
# tensorAdj<-sample_tensor
# impute_nodes(Z_in,tensorAdj,obs_nodes_mat,3)

# OLMF
# OLMF with missing nodes
## The objective function 

lmf_objective<-function(param,tensorAdj,obs_nodes_mat,K){
  n<-dim(tensorAdj)[1]
  L<-dim(tensorAdj)[3]
  
  P<-matrix(param[1:(n*K)],n,K)
  lambda<-lapply(1:L,function(l){return(matrix(param[(n*K+(l-1)*K^2+1):(n*K+l*K^2)],K,K))})
  objloop<- sum(unlist(lapply(1:L,function(l){
    obs_nodes<-obs_nodes_mat[,l]%>%as.logical()
    specobj<-norm(tensorAdj[obs_nodes,obs_nodes,l]-P[obs_nodes,]%*%lambda[[l]]%*%t(P[obs_nodes,]),type="F")^2
    return(specobj)
  })))
  obj=objloop
  return(obj)
}


## The gradient

lmf_grad<-function(param,tensorAdj,obs_nodes_mat,K){
  n<-dim(tensorAdj)[1]
  L<-dim(tensorAdj)[3]
  
  P<-matrix(param[1:(n*K)],n,K)
  lambda<-lapply(1:L,function(l){return(matrix(param[(n*K+(l-1)*K^2+1):(n*K+l*K^2)],K,K))})
  derlist1<-lapply(1:L,function(l){
    obs_nodes<-obs_nodes_mat[,l]%>%as.logical()
    mat<-matrix(0,nrow=n,ncol=K)
    specobj= -2*(tensorAdj[obs_nodes,obs_nodes,l]-P[obs_nodes,]%*%lambda[[l]]%*%t(P[obs_nodes,]))%*%P[obs_nodes,]%*%lambda[[l]]
    mat[obs_nodes,]<-specobj
    return(mat)
  })
  derlist2<-lapply(1:L,function(l){
    obs_nodes<-obs_nodes_mat[,l]%>%as.logical()
    
    specobj= -t(P[obs_nodes,])%*%(tensorAdj[obs_nodes,obs_nodes,l]-P[obs_nodes,]%*%lambda[[l]]%*%t(P[obs_nodes,]))%*%P[obs_nodes,]
    
    return(specobj)
  })
  der1<-Reduce("+",derlist1)
  der2<-unlist(derlist2)
  return(c(as.vector(der1),as.vector(der2)))
}

# The gradient can also be computed on the Stiefel manifold (orthogonal constraint)
lmf_grad2<-function(param,tensorAdj,obs_nodes_mat,K){
  n<-dim(tensorAdj)[1]
  L<-dim(tensorAdj)[3]
  
  P<-matrix(param[1:(n*K)],n,K)
  lambda<-lapply(1:L,function(l){return(matrix(param[(n*K+(l-1)*K^2+1):(n*K+l*K^2)],K,K))})
  derlist1<-lapply(1:L,function(l){
    obs_nodes<-obs_nodes_mat[,l]%>%as.logical()
    mat<-matrix(0,nrow=n,ncol=K)
    specobj= -(diag(n)-P[obs_nodes,]%*%t(P[obs_nodes,]))%*%tensorAdj[obs_nodes,obs_nodes,l]%*%P[obs_nodes,]%*%lambda[[l]]
    mat[obs_nodes,]<-specobj
    return(mat)
  })
  derlist2<-lapply(1:L,function(l){
    obs_nodes<-obs_nodes_mat[,l]%>%as.logical()
    
    specobj= -t(P[obs_nodes,])%*%(tensorAdj[obs_nodes,obs_nodes,l]-P[obs_nodes,]%*%lambda[[l]]%*%t(P[obs_nodes,]))%*%P[obs_nodes,]
    
    return(specobj)
  })
  der1<-Reduce("+",derlist1)
  der2<-unlist(derlist2)
  return(c(as.vector(der1),as.vector(der2)))
}

## BFGS optimization

lmfo_miss<-function(tensorAdj,obs_nodes_mat,K){
  n<-dim(tensorAdj)[1]
  L<-dim(tensorAdj)[3]
  
  # Initialize with mean adjacency
  sum_of_adjmat<-sum_of_adj(tensorAdj)
  spectra<-svds(sum_of_adjmat,K)
  
  P<-spectra$u
  #Same initialization for each lambda_l ?!
  lambda<-lapply(1:L,function(l){return(diag(spectra$d[1:K]))})
  param<-c(as.vector(P),as.vector(unlist(lambda)))
  optimized <-optim(par=param,fn=lmf_objective,gr=lmf_grad,method="BFGS",control=list(reltol=0.0001,maxit=200),
                    tensorAdj=tensorAdj,obs_nodes_mat=obs_nodes_mat,K=K)
  param<-optimized$par
  
  P<-matrix(param[1:(n*K)],n,K)
  lambda<-lapply(1:L,function(l){return(matrix(param[(n*K+(l-1)*K^2+1):(n*K+l*K^2)],K,K))})
  
  specstar<-kmeans(P,K)
  specclus<-specstar$cluster
  return(specclus)
  
}

# Variational EM