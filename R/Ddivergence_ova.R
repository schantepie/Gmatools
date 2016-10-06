Ddivergence_ova <-
function(name,names_pop){
  nb_Gmatrix=length(names_pop)
  name=as.character(name)
  
  ############## Packages requirement
  if (require(mvtnorm) == FALSE)
    stop("mvtnorm not loaded")
  if (require(coda) == FALSE)
    stop("coda not loaded")
  if (require(MCMCglmm) == FALSE)
    stop("coda not loaded")
  
  ############## load first Gmatrix to assess the number of traits and the number of iterations
  model=get(paste(name,1,sep=''))
  # model=model$VCV
  # model=model[,agrep (".animal",colnames(model))]
  nb_trait=sqrt(dim(model)[2])
  nb_iter=dim(model)[1]
  
  ##############Load Gmatrices into an array
  mod=array(,c(nb_iter,nb_trait*nb_trait,nb_Gmatrix)) 
  mod[,,1]=model
  for (i in 2:nb_Gmatrix){ 
    model=get(paste(name,i,sep=''))
    mod[,,i]=model
  }
  siter=sample(nb_iter)
  while (length(which(!siter==1:nb_iter))!= nb_iter) siter=sample(nb_iter)
  sample_mod=mod[siter,,] 
  ################# Function to get distribtion of trait and distance between densities
  Ddivergence<-function (z = NULL, u = NULL, n = 1000) {
    xi <- rmvnorm(n, rep(0, dim(z)[1]),z)
    fx <- dmvnorm(xi, rep(0, dim(z)[1]),z)
    gx <- dmvnorm(xi, rep(0, dim(z)[1]),u)
    mean(sqrt(0.5 * ((fx - gx)^2)/(fx + gx)))
  }
  
  ################ Applying function between the array of Gmatrix and the randomized array
  dist=array(,c(nb_Gmatrix,nb_Gmatrix,nb_iter)) 
  
  
  for (i in 1:nb_iter){
    mo=array(mod[i,,],c(nb_trait,nb_trait,nb_Gmatrix))
    mo_sa=array(sample_mod[i,,],c(nb_trait,nb_trait,nb_Gmatrix))
    for (j in 1:nb_Gmatrix){
      dist[,j,i] =apply(mo,3,function(x,y) {Ddivergence(x,y)},y=mo_sa[,,j]) 
    }    
  }
  
  
  ############## Comparison between intra and intra variation of D
  if(nb_Gmatrix==2){
    Ddiv=mcmc(apply(dist,3, function(x) combn(diag(x), 2,sum,simplify=T)-((x[lower.tri(x)]+t(x)[lower.tri(x)]))))
    Ddist=mcmc(apply(dist,3,function(x) x[lower.tri(x)]))
    difference_summa=cbind(posterior.mode(Ddiv),HPDinterval(Ddiv))
    # Pvalue=length(Ddiv[Ddiv>0])/nb_iter#}
  }else{
    Ddiv=mcmc(t(apply(dist,3, function(x) combn(diag(x), 2,sum,simplify=T)-((x[lower.tri(x)]+t(x)[lower.tri(x)])))))
    Ddist=mcmc(t(apply(dist,3,function(x) x[lower.tri(x)])))
    difference_summa=cbind(posterior.mode(Ddiv),HPDinterval(Ddiv)) 
    # Pvalue=length(Ddiv[Ddiv>0])/nb_iter#}
  }
  Ddistance_summa=cbind(posterior.mode(Ddist),HPDinterval(Ddist))
  names_2by2=apply(do.call(rbind,combn(1:(nb_Gmatrix), 2,simplify=F)),1, paste, collapse="_")
  rownames(difference_summa)=names_2by2
  rownames(Ddistance_summa)=names_2by2
  
  Dmat=matrix("--",nb_Gmatrix,nb_Gmatrix)
  Ddist=matrix(as.character(round(Ddistance_summa,2)),dim(Ddistance_summa)[1],3)
  Dmat[lower.tri(Dmat)]<-c(paste(Ddist[,1]," [",Ddist[,2],":",Ddist[,3],"]",sep=""))
  Dmat=t(Dmat)
  Ddist2=matrix(as.character(round(difference_summa,2)),dim(difference_summa)[1],3)
  Dmat[lower.tri(Dmat)]<-c(paste(Ddist2[,1]," [",Ddist2[,2],":",Ddist2[,3],"]",sep=""))
  Dmat=t(Dmat)
  dimnames(Dmat)=list(names_pop,names_pop)
  results<-list(distance=Ddistance_summa,difference=difference_summa,Results_table=Dmat,Ddist=Ddist)#,pvalue=Pvalue
  return(results)
}
