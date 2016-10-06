krzanowski <-
function(name,names_pop,nb_dimension_subspace=1){
  nb_Gmatrix=length(names_pop)
  names=as.character(name)
  k=nb_dimension_subspace
  
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
  
  if (k>(nb_trait/2))
    stop("k must be < nb_traits/2")
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
  
  
  Krzanowskibytwo<-function (z = NULL, u = NULL,k=k) {
    A <- as.matrix(eigen(z)$vectors[,1:k])
    B <- as.matrix(eigen(u)$vectors[,1:k])
    S=t(A) %*% B %*% t(B) %*% A
    Svectors=eigen(S)$vectors[,1:k] #Eigenvector of S represents the trait combination with the most genetic variation shared among the nine populations.
    Svalues=eigen(S)$values[1:k]     #Eigenvalue can take on a maximum value of 2, which would indicate that the two populations share genetic variance for the trait combination defined by the eigenvector.
    Sum_eigofS <- sum(Svalues)
    angle=acos(Svalues^0.5)/(pi/180)#"/(pi/180)" to get transform angle from radian to degree 
    bi=A%*%Svectors
    MMb=B%*%t(B)%*% bi
    ci=matrix(NA,length(diag(z)),k)
    ci[,1]=(2*(1+sqrt(Svalues[1])))^-0.5* (diag(dim(z)[1]) + ((1/sqrt(Svalues[1]))*(B%*%t(B))))%*% bi[,1]
    for(i in 1:k)  ci[,k]=(2*(1+sqrt(Svalues[k])))^-0.5* (diag(dim(z)[1]) + ((1/sqrt(Svalues[k]))*(B%*%t(B))))%*% bi[,k]
    result<-list(angle=angle,Sum_eigofS=Sum_eigofS,bi=bi,MMb=MMb,ci=ci)
    return(result)
    # return(angle)
  }
  
  ################ Applying function between the array of Gmatrix and the randomized array
  
  dist_angle=array(,c(nb_Gmatrix,nb_Gmatrix,nb_iter,k)) 
  dist_S=array(,c(nb_Gmatrix,nb_Gmatrix,nb_iter,1)) 
  dist_bi=array(,c(nb_trait*nb_Gmatrix,nb_Gmatrix,nb_iter,k)) 
  dist_MMb=array(,c(nb_trait*nb_Gmatrix,nb_Gmatrix,nb_iter,k)) 
  dist_ci=array(,c(nb_trait*nb_Gmatrix,nb_Gmatrix,nb_iter,k)) 
  
   
  for (i in 1:nb_iter){
    mo=array(mod[i,,],c(nb_trait,nb_trait,nb_Gmatrix))
    mo_sa=array(sample_mod[i,,],c(nb_trait,nb_trait,nb_Gmatrix))
    for (j in 1:nb_Gmatrix){
      dist=list()
      dist=apply(mo,3,function(x,y) {Krzanowskibytwo(x,y,k)},y=mo_sa[,,j]) 
      list <- unlist(dist, recursive = FALSE)
      dist_angle[,j,i,] =do.call(rbind,list[grep("angle", names(list))])
      dist_S[,j,i,]=do.call(rbind,list[grep("Sum_eigofS", names(list))])
      dist_bi[,j,i,] =do.call(rbind,list[grep("bi", names(list))])  
      dist_MMb[,j,i,] =do.call(rbind,list[grep("MMb", names(list))])  
      dist_ci[,j,i,] =do.call(rbind,list[grep("ci", names(list))])
    }
  }    
  
  ############## Comparison between intra and intra variation of angle
  
  div=apply(dist_angle,c(3,4), function(x) combn(diag(x), 2,sum,simplify=T)-((x[lower.tri(x)]+t(x)[lower.tri(x)])))
  div_ang=mcmc(t(apply(div,2,c)))
  div_angle_summa=cbind(posterior.mode(mcmc(div_ang)),HPDinterval(mcmc(div_ang)))
  
  angl=apply(dist_angle,c(3,4),function(x) x[lower.tri(x)])
  angl=mcmc(t(apply(angl,2,c)))
  angle_summa=cbind(posterior.mode(angl),HPDinterval(angl))
  
  S=apply(dist_S,c(3,4),function(x) x[lower.tri(x)])
  S=mcmc(t(apply(S,2,c)))
  S_summa=cbind(posterior.mode(S),HPDinterval(S))
  
  
  di_MMb=apply(dist_MMb,c(3,4),function(x) x)
  di_MMb=mcmc(t(apply(di_MMb,2,c)))
  di_MMb_summa=cbind(posterior.mode(di_MMb),HPDinterval(di_MMb))
  di_ci=apply(dist_ci,c(3,4),function(x) x)
  di_ci=mcmc(t(apply(di_ci,2,c)))
  di_ci_summa=cbind(posterior.mode(di_ci),HPDinterval(di_ci))
  di_bi=apply(dist_bi,c(3,4),function(x) x)
  di_bi=mcmc(t(apply(di_bi,2,c)))
  di_bi_summa=cbind(posterior.mode(di_bi),HPDinterval(di_bi))
  
  BM=list()
  for (i in seq(1,(dim(di_bi)[2]),nb_trait)){
    b=di_bi[,i:(i+nb_trait-1)]
    m=di_MMb[,i:(i+nb_trait-1)]
    bm=cbind(b,m)
    bm[which(bm[,1]<0),]=-bm[which(bm[,1]<0),]
    inter=mcmc(bm[,1:nb_trait]-bm[,(nb_trait+1):(nb_trait*2)])  
    BM[[i]]<-cbind(posterior.mode(inter),HPDinterval(inter))
  }
  traits_diff=do.call(rbind,BM)
  
  eig=rep(1:k,each=nb_trait*nb_Gmatrix*nb_Gmatrix)
  trai=rep(1:nb_trait,length.out=length(eig)) 
  ma1=rep(rep(1:nb_Gmatrix,each=nb_trait*nb_Gmatrix),k)
  ma2=rep(rep(1:nb_Gmatrix,each=nb_trait),nb_Gmatrix*k)
  bi_mmb_ci=paste("eig",eig,"_trait",trai,"_",paste(ma1,ma2,sep="-"),sep="")
  
  rownames(di_bi_summa)=bi_mmb_ci
  rownames(di_MMb_summa)=bi_mmb_ci
  rownames(di_ci_summa)=bi_mmb_ci
  rownames(traits_diff)=bi_mmb_ci
  
  up_tri=apply(do.call(rbind,combn(1:(nb_Gmatrix), 2,simplify=F)),1, paste, collapse="_")
  up_tri_k=rep(up_tri,k)
  angle_k=paste(rep("angle",length(up_tri_k)),rep(1:k,each=length(up_tri)),sep="")
  names_2by2_div=paste("divergence",angle_k,up_tri_k,sep="-")
  rownames(div_angle_summa)=names_2by2_div
  
  names_2by2=paste(angle_k,up_tri_k,sep="-")
  rownames(angle_summa)= names_2by2
  rownames(S_summa)= names_2by2
  
  
  Dmat=matrix("--",nb_Gmatrix,nb_Gmatrix)
  Ddist=matrix(as.character(round(div_angle_summa,2)),dim(div_angle_summa)[1],3)
  Dmat[lower.tri(Dmat)]<-c(paste(Ddist[,1]," [",Ddist[,2],":",Ddist[,3],"]",sep=""))
  Dmat=t(Dmat)
  Ddist2=matrix(as.character(round(angle_summa,2)),dim(angle_summa)[1],3)
  Ddist2=c(paste(Ddist2[,1]," [",Ddist2[,2],":",Ddist2[,3],"]",sep=""))
  Dmat[lower.tri(Dmat)]<- Ddist2
  ResMat_Angle=t(Dmat)
  dimnames( ResMat_Angle)=list(names_pop,names_pop)
  
  
  results<-list(
    angle_summa=angle_summa,
    diff_angle_summa=div_angle_summa,
    angle=angl,
    diff_angle=div_ang,
    ResMat_Angle= ResMat_Angle,
    S=S_summa,di_bi_summa=di_bi_summa,
    di_bi_summa=di_bi_summa,
    di_MMb_summa=di_MMb_summa,
    di_ci_summa=di_ci_summa,
    di_bi=di_bi,
    di_MMb=di_MMb,
    traits_diff=traits_diff,
    nb_traits=nb_trait)
  return(results)
}
