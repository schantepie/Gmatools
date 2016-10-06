extract_basic_stats <-
function(name,names_pop,names_traits) {

require(ellipse)
require(MCMCglmm)
require(ggplot2)
  
elli=list()
nb_pop=length(names_pop)
nb_trait=length(names_traits)
combtrait=t(combn(names_traits,2))
coord_elli<-function(x,y,z,k,l){
  MA=matrix(z,k)
  dimnames(MA)=list(names_traits,names_traits)
  MA2=MA[x,x]
  coor=ellipse(MA2,npoints=50)
  tablecoord=as.data.frame(cbind(coor,paste("ite",l,sep=""),paste(x[1],"-",x[2]), paste(names_popi)))
  colnames(tablecoord)=c("x","y","ite","traits","Pop")
  return(tablecoord)
}


#################Table pour le represention des 50 premieres ellipses

allpop_tabl_coor_bycombtrait=list()
for (i in 1:length(names_pop)){
  names_popi=names_pop[i]
  tabl_coor_bycombtrait=list()
  for( d in 1:50){
    Gmat_ij=get(paste(name,i,sep=""))[d,]
    tabl_coor_bycombtrait[[d]]=do.call(rbind,apply(combtrait,1,coord_elli,y=names_popi,z=Gmat_ij,k=nb_trait,l=d))
  }
  allpop_tabl_coor_bycombtrait[[i]]=do.call(rbind,tabl_coor_bycombtrait)
}
Table_plot_ellipse2D_CONFIDENCE=as.data.frame(do.call(rbind,allpop_tabl_coor_bycombtrait))
Table_plot_ellipse2D_CONFIDENCE$x=as.numeric(as.character(Table_plot_ellipse2D_CONFIDENCE$x))
Table_plot_ellipse2D_CONFIDENCE$y=as.numeric(as.character(Table_plot_ellipse2D_CONFIDENCE$y))

#################Table pour le represention du posterieur mode

allpop_tabl_coor_bycombtrait_MODE=list()

for (i in 1:length(names_pop)){
  names_popi=names_pop[i]
  Gmat_ij=posterior.mode(get(paste(name,i,sep="")))
  allpop_tabl_coor_bycombtrait_MODE[[i]]=do.call(rbind,apply(combtrait,1,coord_elli,y=names_popi,z=Gmat_ij,k=nb_trait,l=1))
}
Table_plot_ellipse2D_MODE=as.data.frame(do.call(rbind,allpop_tabl_coor_bycombtrait_MODE))
Table_plot_ellipse2D_MODE$xmode=as.numeric(as.character(Table_plot_ellipse2D_MODE$x))
Table_plot_ellipse2D_MODE$ymode=as.numeric(as.character(Table_plot_ellipse2D_MODE$y))


class(Table_plot_ellipse2D_MODE$ymode)

library(grid)
library(ggplot2)
plot_ellipse = ggplot(data=Table_plot_ellipse2D_CONFIDENCE,aes(x=x,y=y))+
  geom_path(data=Table_plot_ellipse2D_CONFIDENCE,aes(x=x,y=y,group=ite),color="grey")+
  geom_path(data=Table_plot_ellipse2D_MODE, aes(x=xmode,y=ymode),size = 1)+
  facet_grid(Pop~traits)+ 
  scale_x_continuous(name="Additive genetic variance")+
  scale_y_continuous(name="Additive genetic variance")+#"breaks=seq(-2,2,1)
  theme_bw() +
  theme(
    title=element_text(face="bold", size=12),
    panel.background = element_rect(colour="black"),
    panel.margin = unit( 0 ,'mm'),
    axis.ticks=element_line(size=1),
    axis.text=element_text(face="bold", size=12),
    #panel.border=element_rect(size=1,colour="black"),
    axis.title.x = element_text(face="bold", size=12,vjust = -0.4),
    axis.title.y = element_text(face="bold", size=12, angle=90,vjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    legend.position = "none"   
  )

#################""   representation variances et correlations

nb_pop=length(names_pop)
nb_trait=length(names_traits)

mode_CI=array(NA,c(nb_trait*nb_trait+(((nb_trait^2-nb_trait)/2)),3,nb_pop))
for (i in 1:length(names_pop)){
  gmat=get(paste("Gmat_",i,sep=""))
  COVmode_CI=cbind(posterior.mode(gmat),HPDinterval(gmat))
  CORREL_mode_CI=cbind(posterior.mode(posterior.cor(gmat)),HPDinterval(posterior.cor(gmat)))[which(lower.tri(matrix(nrow = nb_trait,ncol = nb_trait))),]
  mode_CI[,,i]= rbind(COVmode_CI,CORREL_mode_CI)
}
names_vacov=apply(expand.grid(names_traits,names_traits),1,function(x) if(x[1]==x[2]) {return (paste("Va(",x[1],")",sep=""))} else {return (paste("Cov(",x[1],",",x[2],")",sep=""))} )
names_cor=apply(t(combn(names_traits,2)),1,function(x) paste("Correl(",x[1],",",x[2],")",sep=""))
dimnames(mode_CI)=list(c(names_vacov,names_cor),c("mode","low","up"),names_pop)

catch_va_correl=c(grep(pattern = "Va",rownames(mode_CI)),grep(pattern = "Correl",rownames(mode_CI)))
return(list(mode_CI=mode_CI,
plot_ellipse=plot_ellipse))
}
