thining_scale <-
function(name=NULL,nb_gmat=NULL,scale=NULL,traits=NULL,div_thin=NULL,mean_std=FALSE,var_std=FALSE, data_ext=NULL){
  require(MCMCglmm)
  
  load(paste(1,name,".Rdata" , sep=""))
  
  if (is.null(data)) data_reduced=data_ext
  
  for (i in 1:length(traits)){
    if(length(agrep (traits[i],colnames(Gmat$VCV)))<1) stop("The names of scaling traits are not in VCV matrix")
  } 
  
  data_reduced=data[,traits]
  # data_reduced=t(apply(data_reduced,1,function(x) {as.numeric(as.character(x))}))
  data_reduced=as.data.frame(data_reduced)
  GMAT_VCV_list=list()
  
  for (i in 1:nb_gmat){
    nomFichier=paste(i,name,".Rdata" , sep="")
    load(nomFichier)
    GMAT_VCV_list[[i]]= assign(paste("Gmat_",i,sep=""),Gmat$VCV)
  }
  GMAT_VCV_list_concat<- do.call(rbind,GMAT_VCV_list)
  

  ############"  subsetting traits
  tri=expand.grid(traits,traits)
  tra=paste("trait",tri[,1],":trait",tri[,2],sep="")
  cut=NULL
  for (i in 1:length(tra)){
    cut=c(cut,grep(tra[i],colnames(GMAT_VCV_list_concat)))
  }
  cut=cut[order(cut)]
  GMAT_VCV_list_concat<-GMAT_VCV_list_concat[,cut]
  ########################" 

  
  gmat_thin=GMAT_VCV_list_concat[,agrep (".animal",colnames(GMAT_VCV_list_concat))]
  if (!is.null(div_thin))  {
    GMAT_VCV_list_concat_thining=dim(GMAT_VCV_list_concat)[1]/div_thin
    GMAT_VCV_list_concat_thin= mcmc(GMAT_VCV_list_concat[seq(1,GMAT_VCV_list_concat_thining*div_thin,div_thin),])
    gmat_thin=GMAT_VCV_list_concat_thin[,grep (".animal",colnames(GMAT_VCV_list_concat))]
  }
  gmat_thin=mcmc(gmat_thin)
  autoc=autocorr(gmat_thin)
  
  grep_pheno<-function(x,y){
    subpart=y[,grep(x,colnames(y))]
    Vp=apply(subpart,1,sum)
    return(Vp)
  }
  Vplist=lapply(tra,function(x,y) grep_pheno(x,y),y=GMAT_VCV_list_concat_thin)
  VP=t(do.call(rbind,Vplist))
  
    
  if (!is.null(scale)){
    data_reduced=t(apply(data_reduced,1,function(x) {as.numeric(x)/(scale)}))
    GMAT_VCV_list_concat=mcmc(t(apply(GMAT_VCV_list_concat[,grep (".animal",colnames(GMAT_VCV_list_concat))],1,function(x) {x/(c(scale%*%t(scale)))})))
    VP=mcmc(t(apply(VP,1,function(x) {x/(c(scale%*%t(scale)))})))
    }

  
  x=as.numeric(data_reduced[,1])
  mean_data=apply(data_reduced,2,function(x)mean(as.numeric(x[!is.na(x)])))
  sd_data=apply(data_reduced,2,function(x)sd(as.numeric(x[!is.na(x)])))
  mean_sd_data=t(rbind(mean_data,sd_data))
  
  
  if (mean_std==TRUE)  {
    tmean_data=t(mean_data)
    std_mean_data=mean_data%*%tmean_data
    colnames(std_mean_data)=colnames(data_reduced)
    rownames(std_mean_data)=colnames(data_reduced)
    std_mean_data=c(std_mean_data)
    gmat_thin=mcmc(t(apply(gmat_thin,1,function(x) {x/std_mean_data})))
    VP_std=mcmc(t(apply(VP,1,function(x) {x/std_mean_data})))
    
  }
  
  if (var_std==TRUE){
    ######one way to standardize by variance is using standard deviation from animal model output  
    tr=paste("trait",traits,sep="")
    tr=paste(tr,":",tr,sep="")
    grep_sd<-function(x,y){
      subpart=y[,grep(x,colnames(y))]
      Standard_deviation_Pheno=sqrt(apply(subpart,1,sum))
    }
    SDlist=lapply(tr,function(x,y) grep_sd(x,y),y=GMAT_VCV_list_concat_thin)
    SD=t(do.call(rbind,  SDlist))
    std_var=t(apply(SD,1,function(x) x%*%t(x)))
    gmat_thin=mcmc(GMAT_VCV_list_concat_thin[,grep(".animal",colnames(GMAT_VCV_list_concat_thin))] / std_var)
    dim(GMAT_VCV_list_concat_thin[,grep(".animal",colnames(GMAT_VCV_list_concat_thin))])
    VP_std=mcmc(VP/ std_var)
    ######the other way to standardize by variance is using standard deviation from row data
    #     sd_data=apply(data_reduced,2,function(x)sd(as.numeric(x[!is.na(x)])))
    #     tsd_data=t(sd_data)
    #     std_var_data=sd_data%*%tsd_data
    #     colnames(std_var_data)=colnames(data_reduced)
    #     rownames(std_var_data)=colnames(data_reduced)
    #     std_var_data=c(std_var_data)
    #     gmat_thin=mcmc(t(apply(gmat_thin,1,function(x) {x/std_var_data})))
    ##############################################"    
    
  }
  return(list(gmat_thin=gmat_thin,autocorr=autoc,mean_sd_data=mean_sd_data,VP=VP,VP_std=VP_std))
}
