krzanowski_vavp <-
function(name,names_pop){
  nb_Gmatrix=length(names_pop)*2
  names=as.character(name)
  cut=length(names_pop)
  ############## Packages requirement
  if (require(mvtnorm) == FALSE)
    stop("mvtnorm not loaded")
  if (require(coda) == FALSE)
    stop("coda not loaded")
  if (require(MCMCglmm) == FALSE)
    stop("coda not loaded")
  
  ############## load first Gmatrix to assess the number of traits and the number of iterations
  model=get(paste(names,1,sep=''))
  # model=model$VCV
  # model=model[,agrep (".animal",colnames(model))]
  nb_trait=sqrt(dim(model)[2])
  nb_iter=dim(model)[1]
  siter=NULL
  ##############Load Gmatrices into an array
  mod=array(,c(nb_iter,nb_trait*nb_trait,nb_Gmatrix)) 
  mod[,,1]=model[1:nb_iter,]
  for (i in 2:nb_Gmatrix){ 
    model=get(paste(names,i,sep=''))
    # model=model$VCV
    # mod[,,i]=model[,agrep (".animal",colnames(model))]
    mod[,,i]=model[1:nb_iter,]
  }
  
  while (length(which(!siter==1:nb_iter))!= nb_iter) siter=sample(nb_iter)
  sample_mod=mod[siter,,]  #randomize array
  
  ################# Function to get distribtion of trait and distance between densities
  angle_first_eigenvector<-function (z = NULL, u = NULL) {
    
    A <- as.matrix(eigen(z)$vectors[,1])
    B <- as.matrix(eigen(u)$vectors[,1])
    S=t(A) %*% B %*% t(B) %*% A
    Svectors=eigen(S)$vectors 
    Svalues=eigen(S)$values    
    angle=acos(Svalues^0.5)/(pi/180)#"/(pi/180)" to get
    return(angle)
    
    }
  
  ################ Applying function between the array of Gmatrix and the randomized array
#   dist=array(,c(nb_Gmatrix,nb_Gmatrix,nb_iter)) 
  dist_P1G1=array(NA,c(nb_iter,cut)) 
  dist_P2G2=array(NA,c(nb_iter,cut)) 
  dist_G1G2=array(NA,c(nb_iter,cut)) 
  dist_P1P2=array(NA,c(nb_iter,cut)) 

  for (i in 1:nb_iter){
    
    G1=array(mod[i,,1:cut],c(nb_trait,nb_trait,cut))
    G2=array(sample_mod[i,,1:cut],c(nb_trait,nb_trait,cut))
    
    P1=array(mod[i,,(cut+1):nb_Gmatrix],c(nb_trait,nb_trait,cut))
    P2=array(sample_mod[i,,(cut+1):nb_Gmatrix],c(nb_trait,nb_trait,cut))

    volstd_G1=array(NA,c(nb_trait,nb_trait,cut))
    volstd_G2=array(NA,c(nb_trait,nb_trait,cut))
    volstd_P1=array(NA,c(nb_trait,nb_trait,cut))
    volstd_P2=array(NA,c(nb_trait,nb_trait,cut))
    
    std_G1=array(NA,c(nb_iter,cut))
    std_G2=array(NA,c(nb_iter,cut))
    std_P1=array(NA,c(nb_iter,cut))
    std_P2=array(NA,c(nb_iter,cut))
    
    
    std_G1[i,]=apply(G1,3, function(x) sum(diag(x)))
    std_G2[i,]=apply(G2,3, function(x) sum(diag(x)))
    std_P1[i,]=apply(P1,3, function(x) sum(diag(x)))
    std_P2[i,]=apply(P2,3, function(x) sum(diag(x)))
    
    for (h in 1:cut){
    volstd_G1[,,h]=G1[,,h]/std_G1[i,h]
    volstd_G2[,,h]=G2[,,h]/std_G2[i,h]
    volstd_P1[,,h]=P1[,,h]/std_P1[i,h]
    volstd_P2[,,h]=P2[,,h]/std_P2[i,h]
    }

    for (j in 1:cut){
      dist_P1G1[i,j]= angle_first_eigenvector(volstd_G1[,,j],volstd_P1[,,j])
      dist_P2G2[i,j]= angle_first_eigenvector(volstd_G2[,,j],volstd_P2[,,j])
      dist_G1G2[i,j]= angle_first_eigenvector(volstd_G1[,,j],volstd_G2[,,j])
      dist_P1P2[i,j]= angle_first_eigenvector(volstd_P1[,,j],volstd_P2[,,j])
    }   
  }
  
  Inter_intra_significant=(dist_G1G2+dist_P1P2)-(dist_P1G1+dist_P2G2)
  colnames(Inter_intra_significant)=names_pop
  summa_Inter_intra_significant=cbind(posterior.mode(mcmc(Inter_intra_significant)),HPDinterval(mcmc(Inter_intra_significant)))
  rownames(summa_Inter_intra_significant)=names_pop
  
  pval=as.matrix(1-apply(Inter_intra_significant,2,function(x,nb_iter) length(x[x<0])/nb_iter,nb_iter=nb_iter))
  colnames(pval)="pvalue"
  
  Distance_angle_PG=dist_P1G1
  colnames(Distance_angle_PG)=names_pop
  summa_Distance_angle_PG=cbind(posterior.mode(mcmc(Distance_angle_PG)),HPDinterval(mcmc(Distance_angle_PG)))
  rownames(summa_Distance_angle_PG)=names_pop
  results<-list(distance=Distance_angle_PG,summary_distance=summa_Distance_angle_PG,summa_difference_angle=summa_Inter_intra_significant,difference_angle=Inter_intra_significant,pvalue=pval)#,pvalue=Pvalue
  return(results)
}
