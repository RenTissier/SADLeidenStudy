# Code by Renaud Tissier

setwd("input_dir_here")

library(coxme)
dataset<-read.table("BasHoogendam_2018_SMRI_GM_dataset",sep="\t",header=T) #importation of the data


################Function needed to compute Kinship matrix from pedigree information#######################
kin_ibd<-function(d2){    
  ufam<-unique(as.vector(d2[,1]))
  nfam<-length(ufam)
  temp<-c(0,0,0,0,0)
  names(temp)<-c("fid", "iid", "pid", "mid", "sex")
  for (i in 1:nfam)
  { 
    print(i)
    pedi <- d2[d2[,1]==ufam[i],1:5]
    names(pedi)<-c("fid", "iid", "pid", "mid", "sex")
    fathers<-unique(pedi[,3])
    mothers<-unique(pedi[,4])
    for (j in fathers)
      if (match(j, pedi[,2], nomatch = 0) == 0 && j!="0"){
        temp<-rbind(temp,data.frame(fid=ufam[i],iid=j,pid=0,mid=0,sex=1,stringsAsFactors=FALSE))}
    for (k in mothers)
      if (match(k, pedi[,2], nomatch = 0) == 0 && k!="0"){
        temp<-rbind(temp,data.frame(fid=ufam[i],iid=k,pid=0,mid=0,sex=2,stringsAsFactors=FALSE))}
    print(pedi)
    print(temp)
    temp<-rbind(temp,pedi)
    
  }
  I <- dim(temp)[1]
  write(t(cbind(temp[2:I,1:5],cbind(rep(0,I-1),rep(0,I-1)))), file="temp.ped", ncol=7)
  write(c("M","MARKER0"), file="temp.dat", ncol=2)    
  write(c("CHR", "MARKER","POSITION"), file="temp.map", ncol=3)            
  write(c("01", "MARKER0",0), file="temp.map", ncol=3, append=T)        
  mer <- system(paste("merlin --ibd -p temp.ped -d temp.dat -m temp.map --bits 250"), intern=F)
  #kin <- read.table("merlin.kin",skip=1)
  #names(kin) <- c("PED","ID1","ID2","POSITION","KINSHIP")
  ibd <- read.table("merlin.ibd",skip=1)
  #names(ibd) <- c("PED","ID1","ID2","MARKER","P0","P1","P2")
  #return(list(kin=kin, ibd=ibd))
  ibd
}

vullphi_merlin<-function(pedi, kini, ibd)
{	
  members<-unique(pedi[,2])
  np<-nrow(pedi)
  #print(members)	
  phi<-matrix(0, np,np)
  for (j in 1:nrow(kini))
  {   
    if (kini[j,2] %in% members && kini[j,3] %in% members){
      k<-which(members %in% kini[j,2])
      #print(k)
      l<-which(members %in% kini[j,3])
      #print(l)	
      phi[k,l]<-kini[j,(ibd+5)]  
      phi[l,k]<-kini[j,(ibd+5)]  
    }
    #print(phi)
  }
  #phin<-phi+t(phi)
  #diag(phin)<-1
  phi
}

rankx<-function(x)
{
  x<-c(x[,1],x[,2])
  ux<-unique(x)
  nux<-length(ux)
  y<-x
  for (i in 1:nux) y[x==ux[i]]<-i
  y<-matrix(y, ncol=2)
  y
}

######################################### Function to compute the kinship matrix for the whole dataset#########
########### 2 parameters
##### dataset : dataset with pedigree information
##### index : vector of index containing IN THE GOOD ORDER the column number of Id family, Id individual, Id father, Id mother and Sex
Kinship.computation<-function(dataset,index){
  dataset2<-dataset[,index]
  family<-unique(dataset2[,1])
  sizeFamilies<-NULL
  varFamList<-vector(mode="list",length=length(family))
  for (i in 1:length(family)){
    dfamily<-dataset2[dataset2[,1]==family[i],]
    sizeFamilies<-c(sizeFamilies,length(dfamily[,1]))
    kinr<-kin_ibd(dfamily)
    kinmat2<-0.5*vullphi_merlin(dfamily,kinr,1) + vullphi_merlin(dfamily,kinr,2)
    varFamList[[i]]<-kinmat2
  }
  kinmat<-matrix(0,nrow=sum(sizeFamilies),ncol=sum(sizeFamilies))
  startingPoint<-c(1)
  for (i in 1:length(sizeFamilies)){
    familySize<-sizeFamilies[i]
    startingIndex<-sum(startingPoint)
    endingIndex<-startingIndex+familySize-1
    kinmat[startingIndex:endingIndex,startingIndex:endingIndex]<-varFamList[[i]]
    startingPoint<-c(startingPoint,familySize)
  }
  return(kinmat)
}

#########################################################################################################

kinmat<-Kinship.computation(dataset[,],c(2,3,5,4,7))  ###computation of the kinship matrix
dimnames(kinmat)<-replicate(2,dataset[,1],simplify=FALSE) ### giving for each column and row of the kinship matrix the proper individual name


automatizedLmekin<-function(dataset,primaryPhenotype,secondaryPhenotype,covariates,kinship){
  for (j in secondaryPhenotype){
    formula<-formula(paste(j," ~ ", 
                           paste(c(primaryPhenotype,covariates), collapse=" + ")," + (1|ppn)"))
    model<-lmekin(formula, data=dataset, varlist=kinmat)
    vcSecondary<-fixef(model)
    results<-cbind(fixef(model),sqrt(diag(vcov(model))),fixef(model)/sqrt(diag(vcov(model))),1-pchisq(fixef(model)^2/sqrt(diag(vcov(model)))^2,1))
    colnames(results)<-c("Value","Std Error","z","p")
    normalTest<-shapiro.test(model$residuals)
    shapiroTest<-c(NA,NA,normalTest$statistic,normalTest$p.value)
    #tablename<-rep(paste(j,"-",primaryPhenotype,".txt",sep=""),(length(covariates)+3))
    outcome<-rep(paste(j,sep=""),(length(covariates)+3))  
    predictor<-rep(paste(primaryPhenotype,sep=""),(length(covariates)+3))
    #print(tablename)
    results<-rbind(results,shapiroTest)
    #print(results)
    #results<-cbind(tablename,results)
    results<-cbind(results,outcome,predictor)
    write.table(results,paste(j,"-",primaryPhenotype,".txt",sep=""),sep="\t")
  }
}


#### Everything with SAD_endophenotype_diagn

# test subcortical volumes non ROIs - (sub)clinical SAD, with interaction term
automatizedLmekin(dataset,primaryPhenotype=c("SAD_endophenotype_diagn"),
                  secondaryPhenotype=c("Rcaud", "Lcaud", "Raccumb", "Laccumb", "Rthal", "Lthal"),
                  covariates=c("Age_demeaned","Sex", "ICV_div1000_demeaned", "SAD_endophenotype_diagn*Age_demeaned", "SAD_endophenotype_diagn*ICV_div1000_demeaned"),kinship=kinmat)

# test cortical surface area non ROIs - (sub)clinical SAD, with interaction term
automatizedLmekin(dataset,primaryPhenotype=c("SAD_endophenotype_diagn"),
                  secondaryPhenotype=c("R_bankssts_surfavg", "L_bankssts_surfavg", "R_cuneus_surfavg", "L_cuneus_surfavg", "R_entorhinal_surfavg", "L_entorhinal_surfavg", 
                                       "R_isthmuscingulate_surfavg", "L_isthmuscingulate_surfavg", "R_lateraloccipital_surfavg", "L_lateraloccipital_surfavg", 
                                       "R_lingual_surfavg", "L_lingual_surfavg", "R_middletemporal_surfavg", "L_middletemporal_surfavg", "R_parahippocampal_surfavg", "L_parahippocampal_surfavg",
                                       "R_paracentral_surfavg", "L_paracentral_surfavg", "R_parsopercularis_surfavg", "L_parsopercularis_surfavg", "R_parsorbitalis_surfavg", "L_parsorbitalis_surfavg",
                                       "R_parstriangularis_surfavg", "L_parstriangularis_surfavg", "R_pericalcarine_surfavg", "L_pericalcarine_surfavg","R_posteriorcingulate_surfavg",
                                       "L_posteriorcingulate_surfavg", "R_frontalpole_surfavg", "L_frontalpole_surfavg" ),
                  covariates=c("Age_demeaned","Sex", "TotalGlobalSurfArea_div1000_demeaned", "SAD_endophenotype_diagn*Age_demeaned", "SAD_endophenotype_diagn*TotalGlobalSurfArea_div1000_demeaned"),kinship=kinmat)

# test cortical thickness non ROIs - (sub)clinical SAD, with interaction term
automatizedLmekin(dataset,primaryPhenotype=c("SAD_endophenotype_diagn"),
                  secondaryPhenotype=c("R_bankssts_thickavg", "L_bankssts_thickavg", "R_cuneus_thickavg", "L_cuneus_thickavg", "R_entorhinal_thickavg", "L_entorhinal_thickavg", 
                                       "R_isthmuscingulate_thickavg", "L_isthmuscingulate_thickavg", "R_lateraloccipital_thickavg", "L_lateraloccipital_thickavg", 
                                       "R_lingual_thickavg", "L_lingual_thickavg", "R_middletemporal_thickavg", "L_middletemporal_thickavg", "R_parahippocampal_thickavg", "L_parahippocampal_thickavg",
                                       "R_paracentral_thickavg", "L_paracentral_thickavg", "R_parsopercularis_thickavg", "L_parsopercularis_thickavg", "R_parsorbitalis_thickavg", "L_parsorbitalis_thickavg",
                                       "R_parstriangularis_thickavg", "L_parstriangularis_thickavg", "R_pericalcarine_thickavg", "L_pericalcarine_thickavg","R_posteriorcingulate_thickavg",
                                       "L_posteriorcingulate_thickavg", "R_frontalpole_thickavg", "L_frontalpole_thickavg" ),
                  covariates=c("Age_demeaned","Sex", "GlobalMeanThickness_demeaned", "SAD_endophenotype_diagn*Age_demeaned", "SAD_endophenotype_diagn*GlobalMeanThickness_demeaned"),kinship=kinmat)

#### Everything with z-score Sa

# test subcortical volumes non ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("Total_Zscore_bothinstruments"),
                  secondaryPhenotype=c("Rcaud", "Lcaud", "Raccumb", "Laccumb", "Rthal", "Lthal"),
                  covariates=c("Age_demeaned","Sex", "ICV_div1000_demeaned"),kinship=kinmat)

# test cortical surface area non ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("Total_Zscore_bothinstruments"),
                  secondaryPhenotype=c("R_bankssts_surfavg", "L_bankssts_surfavg", "R_cuneus_surfavg", "L_cuneus_surfavg", "R_entorhinal_surfavg", "L_entorhinal_surfavg", 
                                       "R_isthmuscingulate_surfavg", "L_isthmuscingulate_surfavg", "R_lateraloccipital_surfavg", "L_lateraloccipital_surfavg", 
                                       "R_lingual_surfavg", "L_lingual_surfavg", "R_middletemporal_surfavg", "L_middletemporal_surfavg", "R_parahippocampal_surfavg", "L_parahippocampal_surfavg",
                                       "R_paracentral_surfavg", "L_paracentral_surfavg", "R_parsopercularis_surfavg", "L_parsopercularis_surfavg", "R_parsorbitalis_surfavg", "L_parsorbitalis_surfavg",
                                       "R_parstriangularis_surfavg", "L_parstriangularis_surfavg", "R_pericalcarine_surfavg", "L_pericalcarine_surfavg","R_posteriorcingulate_surfavg",
                                       "L_posteriorcingulate_surfavg", "R_frontalpole_surfavg", "L_frontalpole_surfavg" ),
                  covariates=c("Age_demeaned","Sex", "TotalGlobalSurfArea_div1000_demeaned"),kinship=kinmat)

# test cortical thickness non ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("Total_Zscore_bothinstruments"),
                  secondaryPhenotype=c("R_bankssts_thickavg", "L_bankssts_thickavg", "R_cuneus_thickavg", "L_cuneus_thickavg", "R_entorhinal_thickavg", "L_entorhinal_thickavg", 
                                       "R_isthmuscingulate_thickavg", "L_isthmuscingulate_thickavg", "R_lateraloccipital_thickavg", "L_lateraloccipital_thickavg", 
                                       "R_lingual_thickavg", "L_lingual_thickavg", "R_middletemporal_thickavg", "L_middletemporal_thickavg", "R_parahippocampal_thickavg", "L_parahippocampal_thickavg",
                                       "R_paracentral_thickavg", "L_paracentral_thickavg", "R_parsopercularis_thickavg", "L_parsopercularis_thickavg", "R_parsorbitalis_thickavg", "L_parsorbitalis_thickavg",
                                       "R_parstriangularis_thickavg", "L_parstriangularis_thickavg", "R_pericalcarine_thickavg", "L_pericalcarine_thickavg","R_posteriorcingulate_thickavg",
                                       "L_posteriorcingulate_thickavg", "R_frontalpole_thickavg", "L_frontalpole_thickavg" ),
                  covariates=c("Age_demeaned","Sex", "GlobalMeanThickness_demeaned"),kinship=kinmat)

#### Everything with FNE

# test subcortical volumes non ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("FNE"),
                  secondaryPhenotype=c("Rcaud", "Lcaud", "Raccumb", "Laccumb", "Rthal", "Lthal"),
                  covariates=c("Age_demeaned","Sex", "ICV_div1000_demeaned"),kinship=kinmat)

# test cortical surface area non ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("FNE"),
                  secondaryPhenotype=c("R_bankssts_surfavg", "L_bankssts_surfavg", "R_cuneus_surfavg", "L_cuneus_surfavg", "R_entorhinal_surfavg", "L_entorhinal_surfavg", 
                                       "R_isthmuscingulate_surfavg", "L_isthmuscingulate_surfavg", "R_lateraloccipital_surfavg", "L_lateraloccipital_surfavg", 
                                       "R_lingual_surfavg", "L_lingual_surfavg", "R_middletemporal_surfavg", "L_middletemporal_surfavg", "R_parahippocampal_surfavg", "L_parahippocampal_surfavg",
                                       "R_paracentral_surfavg", "L_paracentral_surfavg", "R_parsopercularis_surfavg", "L_parsopercularis_surfavg", "R_parsorbitalis_surfavg", "L_parsorbitalis_surfavg",
                                       "R_parstriangularis_surfavg", "L_parstriangularis_surfavg", "R_pericalcarine_surfavg", "L_pericalcarine_surfavg","R_posteriorcingulate_surfavg",
                                       "L_posteriorcingulate_surfavg", "R_frontalpole_surfavg", "L_frontalpole_surfavg" ),
                  covariates=c("Age_demeaned","Sex", "TotalGlobalSurfArea_div1000_demeaned"),kinship=kinmat)

# test cortical thickness non ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("FNE"),
                  secondaryPhenotype=c("R_bankssts_thickavg", "L_bankssts_thickavg", "R_cuneus_thickavg", "L_cuneus_thickavg", "R_entorhinal_thickavg", "L_entorhinal_thickavg", 
                                       "R_isthmuscingulate_thickavg", "L_isthmuscingulate_thickavg", "R_lateraloccipital_thickavg", "L_lateraloccipital_thickavg", 
                                       "R_lingual_thickavg", "L_lingual_thickavg", "R_middletemporal_thickavg", "L_middletemporal_thickavg", "R_parahippocampal_thickavg", "L_parahippocampal_thickavg",
                                       "R_paracentral_thickavg", "L_paracentral_thickavg", "R_parsopercularis_thickavg", "L_parsopercularis_thickavg", "R_parsorbitalis_thickavg", "L_parsorbitalis_thickavg",
                                       "R_parstriangularis_thickavg", "L_parstriangularis_thickavg", "R_pericalcarine_thickavg", "L_pericalcarine_thickavg","R_posteriorcingulate_thickavg",
                                       "L_posteriorcingulate_thickavg", "R_frontalpole_thickavg", "L_frontalpole_thickavg"),
                  covariates=c("Age_demeaned","Sex", "GlobalMeanThickness_demeaned"),kinship=kinmat)

########################################################################################################################
