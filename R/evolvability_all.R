evolvability_all <-function(name,names_pop,cpus=1,simu=FALSE,analytic=TRUE){
  
  require(MCMCglmm)
  require(snow)
  nb_pop=length(names_pop)

  
  results=list()
  summa_EVOL_MEASURES_beta=array(NA,c(nb_pop,3,5))
  summa_EVOL_MEASURES_analytic=array(NA,c(nb_pop,3,5))
  EVOL_MEASURES_analytic=array(NA,c(1000,5,nb_pop))
  
  
for(i in 1:nb_pop){
  gmat=get(paste(name,i,sep=""))
  numtr=sqrt(dim(gmat)[2])
############# estimation with beta
if(simu==TRUE) {

  #beta generation
trial=1000
beta = matrix(rnorm(trial*numtr), nrow = numtr)
beta = t(beta)/sqrt(colSums(beta^2)) 

Cond_evol_resp_auto_inte<-function(x) {
  G=matrix(x,numtr,numtr)
  Ginv=solve(G)
  G2=G%*%G
  eProj=apply(beta,1,function(x) {t(x)%*%G%*%x})
  eC=apply(beta,1,function(x) {1/(t(x)%*%Ginv%*%x)})
  eR=apply(beta,1,function(x) {sqrt((t(x)%*%G2%*%x))})
  Auto = eC/eProj
  Inte = 1-Auto
  mean_eU=mean(eProj/numtr)
  mean_eC=mean(eC)
  mean_eR=mean(eR)
  mean_Auto=mean(Auto)
  mean_Inte=mean(Inte)
  return(list(mean_eU=mean_eU,
              mean_eC=mean_eC,
              mean_eR=mean_eR,
              mean_Auto=mean_Auto,
              mean_Inte=mean_Inte))
}

  clus <- makeCluster(cpus)
  clusterExport(clus, list("beta","numtr"), envir=environment())
  Evol_measure=parApply(clus,gmat,1,Cond_evol_resp_auto_inte)
  stopCluster(clus)
  EVOL_MEASURE=as.matrix(t(simplify2array(Evol_measure,higher = FALSE)))
  EVOL_MEASURE=mcmc( t(apply(EVOL_MEASURE,1,function(x) simplify2array(x, higher = TRUE))))
  MEASURE=cbind(posterior.mode(EVOL_MEASURE),HPDinterval(EVOL_MEASURE))
  summa_EVOL_MEASURES_beta[i,,]=t(MEASURE)
  dimnames(summa_EVOL_MEASURES_beta)=list(names_pop,c("mode","lower","upper"),c("Unconditional_Evolvability","Contional_Evolvability","Respondability","Autonomy","Integration"))
  results[["summa_EVOL_MEASURES_beta"]]=summa_EVOL_MEASURES_beta
}
  
############# analytic estimation
if( analytic==TRUE){  
    
    Cond_evol_analytic<-function(x) {
      
      G=matrix(x,numtr,numtr)
      eig=eigen(G)
      H=1/(mean(1/eig$values))
      Iinv=var(1/eig$values)/(mean(1/eig$values)^2)
      I2=var(eig$values^2)/((mean(eig$values^2))^2)
      I=var(eig$values)/((mean(eig$values))^2)
      k= dim(G)[1]
      
      mean_eU=mean(eig$values)
      mean_eC=H*(1+((2*Iinv)/(k+2)))
      mean_eR = sqrt(mean(eig$values^2))*(1-(I2/(4*(k+2))))
      mean_Auto = (H/mean(eig$values))*(1+2*(I+Iinv-1+(H/mean(eig$values))+2*I*Iinv/(k+2))/(k+2))
      mean_Inte = 1-mean_Auto

      return(list(mean_eU=mean_eU,
                  mean_eC=mean_eC,
                  mean_eR=mean_eR,
                  mean_Auto=mean_Auto,
                  mean_Inte=mean_Inte))
    }
    
    clus <- makeCluster(cpus)
    
    clusterExport(clus, list("numtr"), envir=environment())
    Evol_measure=parApply(clus,gmat,1, Cond_evol_analytic)
    stopCluster(clus)
    EVOL_MEASURE=as.matrix(t(simplify2array(Evol_measure,higher = FALSE)))
    EVOL_MEASURE=mcmc( t(apply(EVOL_MEASURE,1,function(x) simplify2array(x, higher = TRUE))))
    MEASURE=cbind(posterior.mode(EVOL_MEASURE),HPDinterval(EVOL_MEASURE))
    summa_EVOL_MEASURES_analytic[i,,]=t(MEASURE)
    EVOL_MEASURES_analytic[,,i]=EVOL_MEASURE
    dimnames( summa_EVOL_MEASURES_analytic)=list(names_pop,c("mode","lower","upper"),c("Unconditional_Evolvability","Conditional_Evolvability","Respondability","Autonomy","Integration"))
    results[["summa_EVOL_MEASURES_analytic"]]=list(summa_EVOL_MEASURES_analytic=summa_EVOL_MEASURES_analytic,EVOL_MEASURES_analytic=EVOL_MEASURES_analytic)
    }
  }
return (results)
}




