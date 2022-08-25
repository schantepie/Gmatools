shape <-
function(name,names_pop){
  nb_Gmatrix=length(names_pop)
  names=as.character(name)
  
  
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
  
  ##############Load Gmatrices into an array
  mod=array(,c(nb_iter,nb_trait*nb_trait,nb_Gmatrix)) 
  mod[,,1]=model
  for (i in 2:nb_Gmatrix){ 
    model=get(paste(names,i,sep=''))
    # model=model$VCV
    # mod[,,i]=model[,agrep (".animal",colnames(model))]
    mod[,,i]=model
  }
  
  siter=sample(nb_iter)
  while (length(which(!siter==1:nb_iter))!= nb_iter) siter=sample(nb_iter)
  sample_mod=mod[siter,,]  #randomize array
  
  ################# Function to get Krzanowski angles between subspaces
  
  shape_matrix<-function (z = NULL, u = NULL) {
    
    A <-eigen(z)
    B <-eigen(u)
    
    Vol_traceA=sum(diag(z))
    Vol_traceB=sum(diag(u))
    dif_vol_trace=Vol_traceA-Vol_traceB
    
    Vol_piA=(1/2)*pi^2*prod(sqrt(A$values))
    Vol_piB=(1/2)*pi^2*prod(sqrt(B$values))
    dif_vol_pi=Vol_piA-Vol_piB
    
    EccentricityA=A$values[1]/A$values[2]
    EccentricityB=B$values[1]/B$values[2]
    dif_Eccentricity=EccentricityA-EccentricityB
    
    Var_maxA=abs(A$values[1])/sum(abs(A$values))
    Var_maxB=abs(B$values[1])/sum(abs(B$values))
    dif_Var_max=Var_maxA-Var_maxB
    #Evol_random=
    result<-list(
      dif_vol_trace=dif_vol_trace,
      dif_vol_pi=dif_vol_pi,
      dif_Eccentricity=dif_Eccentricity,
      dif_Var_max=dif_Var_max)
    return(result)
    # return(angle)
  }
  
  ################ Applying function between the array of Gmatrix and the randomized array
  
  dif_vol_trace=array(,c(nb_Gmatrix,nb_Gmatrix,nb_iter)) 
  dif_vol_pi=array(,c(nb_Gmatrix,nb_Gmatrix,nb_iter)) 
  dif_Eccentricity=array(,c(nb_Gmatrix,nb_Gmatrix,nb_iter)) 
  dif_Var_max=array(,c(nb_Gmatrix,nb_Gmatrix,nb_iter)) 
  
  for (i in 1:nb_iter){
    mo=array(mod[i,,],c(nb_trait,nb_trait,nb_Gmatrix))
    mo_sa=array(sample_mod[i,,],c(nb_trait,nb_trait,nb_Gmatrix))
    for (j in 1:nb_Gmatrix){
      dist=list()
      dist=apply(mo,3,function(x,y) {shape_matrix(x,y)},y=mo_sa[,,j]) 
      list <- unlist(dist, recursive = FALSE)
      dif_vol_trace[,j,i] =do.call(rbind,list[grep("dif_vol_trace", names(list))])
      dif_vol_pi[,j,i] =do.call(rbind,list[grep("dif_vol_pi", names(list))])
      dif_Eccentricity[,j,i] =do.call(rbind,list[grep("dif_Eccentricity", names(list))])
      dif_Var_max[,j,i] =do.call(rbind,list[grep("dif_Var_max", names(list))])
    }
  }    
  
  ############## Comparison between intra and intra variation of angle
  
  vol_trace=apply(dif_vol_trace,3, function(x) combn(diag(x), 2,sum,simplify=T)-((x[lower.tri(x)]+t(x)[lower.tri(x)])))
  vol_trace=mcmc(t(apply(vol_trace,2,c)))
  vol_trace_summa=cbind(posterior.mode(mcmc(vol_trace)),HPDinterval(mcmc(vol_trace)))
  # Pvalue_vol_trace=1-apply(vol_trace,2,function(x) length(which(x<0)))/nb_iter
  diffvoltrace=apply(dif_vol_trace,3,function(x) x[lower.tri(x)])
  diffvoltrace=mcmc(t(apply(diffvoltrace,2,c)))
  diffvoltrace_summa=cbind(posterior.mode(diffvoltrace),HPDinterval(diffvoltrace))
  
  
  vol_pi=apply(dif_vol_pi,3, function(x) combn(diag(x), 2,sum,simplify=T)-((x[lower.tri(x)]+t(x)[lower.tri(x)])))
  vol_pi=mcmc(t(apply(vol_pi,2,c)))
  vol_pi_summa=cbind(posterior.mode(mcmc(vol_pi)),HPDinterval(mcmc(vol_pi)))
  # Pvalue_vol_pi=1-apply(vol_pi,2,function(x) length(which(x<0)))/nb_iter
  diffvol_pi=apply(dif_vol_pi,3,function(x) x[lower.tri(x)])
  diffvol_pi=mcmc(t(apply(diffvol_pi,2,c)))
  diffvol_pi_summa=cbind(posterior.mode(diffvol_pi),HPDinterval(diffvol_pi))
  
  
  Eccentricity=apply(dif_Eccentricity,3, function(x) combn(diag(x), 2,sum,simplify=T)-((x[lower.tri(x)]+t(x)[lower.tri(x)])))
  Eccentricity=mcmc(t(apply(Eccentricity,2,c)))
  Eccentricity_summa=cbind(posterior.mode(mcmc(Eccentricity)),HPDinterval(mcmc(Eccentricity)))
  # Pvalue_Eccentricity=1-apply(Eccentricity,2,function(x) length(which(x<0)))/nb_iter
  diffEccentricity=apply(dif_Eccentricity,3,function(x) x[lower.tri(x)])
  diffEccentricity=mcmc(t(apply(diffEccentricity,2,c)))
  diffEccentricity_summa=cbind(posterior.mode(diffEccentricity),HPDinterval(diffEccentricity))
  
  Var_gmax=apply(dif_Var_max,3, function(x) combn(diag(x), 2,sum,simplify=T)-((x[lower.tri(x)]+t(x)[lower.tri(x)])))
  Var_gmax=mcmc(t(apply(Var_gmax,2,c)))
  Var_gmax_summa=cbind(posterior.mode(mcmc(Var_gmax)),HPDinterval(mcmc(Var_gmax)))
  # Pvalue_Var_gmax=1-apply(Var_gmax,2,function(x) length(which(x<0)))/nb_iter
  diffVar_gmax=apply(dif_Var_max,3,function(x) x[lower.tri(x)])
  diffVar_gmax=mcmc(t(apply(diffVar_gmax,2,c)))
  diffVar_gmax_summa=cbind(posterior.mode(diffVar_gmax),HPDinterval(diffVar_gmax))
  
  
  names_2by2=apply(do.call(rbind,combn(1:(nb_Gmatrix), 2,simplify=F)),1, paste, collapse="_")
  rownames(diffvoltrace_summa)= names_2by2
  rownames(diffvol_pi_summa)= names_2by2
  rownames(diffEccentricity_summa)= names_2by2
  rownames(diffVar_gmax_summa)= names_2by2
  rownames(vol_trace_summa)= names_2by2
  rownames(vol_pi_summa)= names_2by2
  rownames(Eccentricity_summa)= names_2by2
  rownames(Var_gmax_summa)=names_2by2
  
  
  Dmat=matrix("--",nb_Gmatrix,nb_Gmatrix)
  Ddist=matrix(as.character(round(diffvoltrace_summa,2)),dim(diffvoltrace_summa)[1],3)
  Dmat[lower.tri(Dmat)]<-c(paste(Ddist[,1]," [",Ddist[,2],":",Ddist[,3],"]",sep=""))
#   Dmat=t(Dmat)
#   Ddist2=matrix(as.character(round(vol_trace_summa,2)),dim(vol_trace_summa)[1],3)
#   Ddist2=c(paste(Ddist2[,1]," [",Ddist2[,2],":",Ddist2[,3],"]",sep=""))
#   Dmat[lower.tri(Dmat)]<- Ddist2
  ResMat_volumes=Dmat
  dimnames(ResMat_volumes)=list(names_pop,names_pop)
  
  Dmat=matrix("--",nb_Gmatrix,nb_Gmatrix)
  Ddist=matrix(as.character(round(diffEccentricity_summa,2)),dim(diffEccentricity_summa)[1],3)
  Dmat[lower.tri(Dmat)]<-c(paste(Ddist[,1]," [",Ddist[,2],":",Ddist[,3],"]",sep=""))
#   Dmat=t(Dmat)
#   Ddist2=matrix(as.character(round(Eccentricity_summa,2)),dim(Eccentricity_summa)[1],3)
#   Ddist2=c(paste(Ddist2[,1]," [",Ddist2[,2],":",Ddist2[,3],"]",sep=""))
#   Dmat[lower.tri(Dmat)]<- Ddist2
 ResMat_eccentricity=Dmat
  dimnames(ResMat_eccentricity)=list(names_pop,names_pop)
  
  Dmat=matrix("--",nb_Gmatrix,nb_Gmatrix)
  Ddist=matrix(as.character(round(diffVar_gmax_summa,2)),dim(diffVar_gmax_summa)[1],3)
  Dmat[lower.tri(Dmat)]<-c(paste(Ddist[,1]," [",Ddist[,2],":",Ddist[,3],"]",sep=""))
#   Dmat=t(Dmat)
#   Ddist2=matrix(as.character(round(Var_gmax_summa,2)),dim(Var_gmax_summa)[1],3)
#   Ddist2=c(paste(Ddist2[,1]," [",Ddist2[,2],":",Ddist2[,3],"]",sep=""))
#   Dmat[lower.tri(Dmat)]<- Ddist2
ResMat_Vargmax=(Dmat)
  dimnames(ResMat_Vargmax)=list(names_pop,names_pop)

  
  results<-list(
#     diffvoltrace_summa=diffvoltrace_summa,
#     diffvol_pi_summa=diffvol_pi_summa,
#     diffEccentricity_summa=diffEccentricity_summa,
#     diffVar_gmax_summa=diffVar_gmax_summa,
#     ConfInterval_vol_trace= vol_trace_summa,
#     ConfInterval_vol_pi= vol_pi_summa,
#     ConfInterval_Eccentricity= Eccentricity_summa,
#     ConfIntervalVar_gmax = Var_gmax_summa,
    ResMat_volumes=ResMat_volumes,
    ResMat_eccentricity=ResMat_eccentricity,
    ResMat_Vargmax=ResMat_Vargmax
    )
  return(results)
}
