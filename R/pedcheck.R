pedcheck <-
function(Ped_check,names_remove=NULL,remove_punk=F,delete_space=F,data_id=NULL){
  
  require(pedigree)
  require(MasterBayes)
  
  Ped_check[,1]=as.character(Ped_check[,1])
  Ped_check[,2]=as.character(Ped_check[,2])
  Ped_check[,3]=as.character(Ped_check[,3])
  
  #-----------------resolve loop of individuals (noted chick and parent in the same line)----------------
  options(warn = 0)
  ind_loop=pedigree::orderPed(Ped_check)#indiv with -1 has a problem
  loop_ind="no loop with individuals"
  if(length(which(ind_loop==-1))>0){
    loop_ind=Ped_check[which(ind_loop==-1),]
    Ped_check=Ped_check[-which(ind_loop==-1),]
    
  }
  
  #----------------delete space in names-------------------------------------------------------------------
  
  if (delete_space==T){
    Ped_check[,1]=gsub(" ", "", Ped_check[,1], fixed = TRUE)
    Ped_check[,2]=gsub(" ", "", Ped_check[,2], fixed = TRUE)
    Ped_check[,3]=gsub(" ", "", Ped_check[,3], fixed = TRUE)
  }

  
  #----------------"delete row with "XXXX given exotic name" inside---------------------
  
  if (!is.null(names_remove)) {
    #one or several text pattern can be used at the same time
    if(length(which(Ped_check [,1]%in%names_remove))>0) Ped_check=Ped_check[-which(Ped_check [,1]%in%names_remove),]
    if(length(which(Ped_check [,2]%in%names_remove))>0) Ped_check[which(Ped_check [,2]%in%names_remove), 2]=NA
    if(length(which(Ped_check [,3]%in%names_remove))>0) Ped_check[which(Ped_check [,3]%in%names_remove), 3]=NA
  }
  
  #---------------- find exotic names remaining
  
  exotic_names_id=NULL
  Ped_check_line=c(as.matrix(Ped_check))
  exotic_names=unique(Ped_check_line[grep("[[:punct:]]", Ped_check_line)])
  
  if (remove_punk==F){
    exotic_names_id="no names with punctuation were detected"
    if (length(exotic_names)>0)      exotic_names_id=exotic_names 
  }
  
  if (remove_punk==T){
    exotic_names_id="no names with punctuation were detected"
    if (length(exotic_names)>0){
      exotic_names_id=exotic_names
      if(length(which(Ped_check [,1]%in%exotic_names_id))>0) Ped_check=Ped_check[-which(Ped_check [,1]%in%exotic_names_id),]
      if(length(which(Ped_check [,2]%in%exotic_names_id))>0) Ped_check[which(Ped_check [,2]%in%exotic_names_id), 2]=NA
      if(length(which(Ped_check [,3]%in%exotic_names_id))>0) Ped_check[which(Ped_check [,3]%in%exotic_names_id), 3]=NA
    }
  }
  
  #-------------------"resolve loop in family (part in pogress) but break if problem-------------------------------
  
loop_family="no family loop"
  
fam_loop<-function(id.parent,id){
  depth=id
  count=0
  while (!all(is.na(id.parent))){
    count= count+1
    id.pass=match(id.parent,id,nomatch= NA)
    depth=cbind(depth,id.parent[id.pass])
    id.parent=depth[,(count+1)]
    if(count==50) {
      cat("familly loop with parents")
      depth=depth[!is.na(depth[,50]),]
      break
    }
  }
  results=list(depth=depth) 
  return(results)
}

mom_loop=fam_loop(Ped_check[,3],Ped_check[,1])
dad_loop=fam_loop(Ped_check[,2],Ped_check[,1])
loop=list(mom_loop,dad_loop)
depth_ped=max(dim(mom_loop$depth)[2],dim(dad_loop$depth)[2])-1
if(depth_ped==50)loop_family=loop

#-------------------get and delete birds noted as dam and sire-------------------------------------------------
  
  pb_damsire= cbind(1:length(match(Ped_check[,2], Ped_check[,3])),match(Ped_check[,2], Ped_check[,3],nomatch=-666,incomparables=NA))
  pb_damsire= pb_damsire[which(pb_damsire[,2]>1),]#if empty no problem
  
  ind_damsire="no dam noted also as a sire"
  if (length(c(pb_damsire))>0){
    pb_damsire= pb_damsire[complete.cases(pb_damsire[,1]),]
    pb_damsire=cbind(pb_damsire,Ped_check[pb_damsire[,1],2])#get id name with probkem on the female column
    ind_damsire=unique(pb_damsire[,3]) #give unique name with problem
    #delete row with problem
    if(length(which(Ped_check [,1]%in%ind_damsire))>0) Ped_check=Ped_check[-which(Ped_check [,1]%in%ind_damsire),]
    if(length(which(Ped_check [,2]%in%ind_damsire))>0) Ped_check=Ped_check[-which(Ped_check [,2]%in%ind_damsire),]
    if(length(which(Ped_check [,3]%in%ind_damsire))>0) Ped_check=Ped_check[-which(Ped_check [,3]%in%ind_damsire),]
  }
  
  
  #------------------------------find duplicated row in pedigree and add founders------------------------------
  
  Ped_check_allfound= Ped_check
  colnames(Ped_check_allfound)=c("id","dam","sire")
  id_female_founder=as.character(Ped_check[,2])
  id_female_founder=id_female_founder[!is.na(id_female_founder)]
  id_female_founder=id_female_founder[which(match(id_female_founder,as.character(Ped_check[,1]),nomatch =-666)<0)]
  if(length(id_female_founder)>0){
  id_female_founder=cbind(id_female_founder,NA,NA)
  colnames(id_female_founder)=c("id","dam","sire")
  Ped_check_allfound=rbind(Ped_check_allfound,id_female_founder)
  }
  
  id_male_founder=as.character(Ped_check[,3])
  id_male_founder=id_male_founder[!is.na(id_male_founder)]
  id_male_founder=id_male_founder[which(match(id_male_founder,as.character(Ped_check[,1]),nomatch =-666)<0)]
  if(length(id_male_founder)>0){
  id_male_founder=cbind(id_male_founder,NA,NA)
  colnames(id_male_founder)=c("id","dam","sire")
  Ped_check_allfound=rbind(Ped_check_allfound,id_male_founder)
  }
  
  Ped_check_allfound=Ped_check_allfound[!duplicated(Ped_check_allfound),]
  
  
  #----check for similarity between data-id and ped-id (if data_id is given individuals missing are added)----------
  
  if (!is.null(data_id)) {
    similar_ped_data_sire=cbind(unique(as.character(data_id[which(match(as.character(data_id),as.character(Ped_check_allfound[,1]),nomatch =-666)<0)])),"NA","NA")
   # if (all(is.na(similar_ped_data_sire))){ 
    colnames(similar_ped_data_sire)=c("id","dam","sire")
    Ped_check_allfound=rbind(Ped_check_allfound,similar_ped_data_sire)
  #  }
  }
  
  #------------------------------------------check for multiparentality------------------------------------------
  
  multiparent="no multiparentality"
  
  if(length(which(duplicated(Ped_check_allfound[,1])))>0){
    multiparent=Ped_check_allfound[duplicated(Ped_check_allfound[,1]),]
  }
  
  #-----------------------------------------final ordering using MasterBayes--------------------------------------
  
  Ped_check_allfound_ord=MasterBayes::orderPed(Ped_check_allfound)
  
  results<-list(
    ped_checked=Ped_check_allfound_ord,
    loop_ind=loop_ind,
    loop_family=loop_family,
    depth_ped=depth_ped,
    ancestors=loop,
    ind_damsire=ind_damsire,
    exotic_names_id=exotic_names_id,
    multiparentality=multiparent)
  return(results)
}
