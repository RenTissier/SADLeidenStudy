# Code by Renaud Tissier

setwd("input_dir_here")

library(coxme)
dataset<-read.table("BasHoogendam_2018_SMRI_GM_dataset.txt",sep="\t",header=T) #importation of the data


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
# test basic imaging phenotypes 
automatizedLmekin(dataset,primaryPhenotype=c("SAD_endophenotype_diagn"),secondaryPhenotype=c("TotalGlobalSurfArea_div1000_demeaned", "ICV_div1000_demeaned", "GlobalMeanThickness_demeaned"),covariates=c("Age_demeaned","Sex"),kinship=kinmat)

# test subcortical volumes ROIs - (sub)clinical SAD, with interaction term
automatizedLmekin(dataset,primaryPhenotype=c("SAD_endophenotype_diagn"),
                  secondaryPhenotype=c("Ramyg","Lamyg","Rhippo", "Lhippo", "Rpal", "Lpal", "Rput", "Lput"),
                  covariates=c("DEP_Z","Age_demeaned","Sex", "ICV_div1000_demeaned", "SAD_endophenotype_diagn*Age_demeaned", "SAD_endophenotype_diagn*ICV_div1000_demeaned"),kinship=kinmat)

# test cortical surface area ROIs - (sub)clinical SAD, with interaction term
automatizedLmekin(dataset,primaryPhenotype=c("SAD_endophenotype_diagn"),
                  secondaryPhenotype=c("R_superiorfrontal_surfavg","L_superiorfrontal_surfavg","R_caudalmiddlefrontal_surfavg", "L_caudalmiddlefrontal_surfavg", "R_rostralmiddlefrontal_surfavg", "L_rostralmiddlefrontal_surfavg", "R_lateralorbitofrontal_surfavg", "L_lateralorbitofrontal_surfavg",
                                       "R_medialorbitofrontal_surfavg", "L_medialorbitofrontal_surfavg", "R_precentral_surfavg", "L_precentral_surfavg",
                                       "R_caudalanteriorcingulate_surfavg", "L_caudalanteriorcingulate_surfavg", "R_rostralanteriorcingulate_surfavg", "L_rostralanteriorcingulate_surfavg",
                                       "R_rostralanteriorcingulate_surfavg", "L_rostralanteriorcingulate_surfavg", "R_insula_surfavg", "L_insula_surfavg",
                                       "R_superiorparietal_surfavg", "L_superiorparietal_surfavg","R_inferiorparietal_surfavg", "L_inferiorparietal_surfavg",
                                       "R_precuneus_surfavg", "L_precuneus_surfavg", "R_supramarginal_surfavg", "L_supramarginal_surfavg",
                                       "R_postcentral_surfavg", "L_postcentral_surfavg", "R_temporalpole_surfavg", "L_temporalpole_surfavg", 
                                       "R_inferiortemporal_surfavg", "L_inferiortemporal_surfavg", "R_superiortemporal_surfavg", "L_superiortemporal_surfavg",
                                       "R_fusiform_surfavg", "L_fusiform_surfavg", "R_transversetemporal_surfavg", "L_transversetemporal_surfavg"),
                  covariates=c("DEP_Z","Age_demeaned","Sex", "TotalGlobalSurfArea_div1000_demeaned", "SAD_endophenotype_diagn*Age_demeaned", "SAD_endophenotype_diagn*TotalGlobalSurfArea_div1000_demeaned"),kinship=kinmat)

# test cortical thickness ROIs - (sub)clinical SAD, with interaction term
automatizedLmekin(dataset,primaryPhenotype=c("SAD_endophenotype_diagn"),
                  secondaryPhenotype=c("R_superiorfrontal_thickavg","L_superiorfrontal_thickavg","R_caudalmiddlefrontal_thickavg", "L_caudalmiddlefrontal_thickavg", "R_rostralmiddlefrontal_thickavg", "L_rostralmiddlefrontal_thickavg", "R_lateralorbitofrontal_thickavg", "L_lateralorbitofrontal_thickavg",
                                       "R_medialorbitofrontal_thickavg", "L_medialorbitofrontal_thickavg", "R_precentral_thickavg", "L_precentral_thickavg",
                                       "R_caudalanteriorcingulate_thickavg", "L_caudalanteriorcingulate_thickavg", "R_rostralanteriorcingulate_thickavg", "L_rostralanteriorcingulate_thickavg",
                                       "R_rostralanteriorcingulate_thickavg", "L_rostralanteriorcingulate_thickavg", "R_insula_thickavg", "L_insula_thickavg",
                                       "R_superiorparietal_thickavg", "L_superiorparietal_thickavg","R_inferiorparietal_thickavg", "L_inferiorparietal_thickavg",
                                       "R_precuneus_thickavg", "L_precuneus_thickavg", "R_supramarginal_thickavg", "L_supramarginal_thickavg",
                                       "R_postcentral_thickavg", "L_postcentral_thickavg", "R_temporalpole_thickavg", "L_temporalpole_thickavg", 
                                       "R_inferiortemporal_thickavg", "L_inferiortemporal_thickavg", "R_superiortemporal_thickavg", "L_superiortemporal_thickavg",
                                       "R_fusiform_thickavg", "L_fusiform_thickavg", "R_transversetemporal_thickavg", "L_transversetemporal_thickavg"),
                  covariates=c("DEP_Z","Age_demeaned","Sex", "GlobalMeanThickness_demeaned", "SAD_endophenotype_diagn*Age_demeaned", "SAD_endophenotype_diagn*GlobalMeanThickness_demeaned"),kinship=kinmat)

#### Everything with z-score Sa
# test basic imaging phenotypes 
automatizedLmekin(dataset,primaryPhenotype=c("Total_Zscore_bothinstruments"),secondaryPhenotype=c("TotalGlobalSurfArea_div1000_demeaned", "ICV_div1000_demeaned", "GlobalMeanThickness_demeaned"),covariates=c("Age_demeaned","Sex"),kinship=kinmat)

# test subcortical volumes ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("Total_Zscore_bothinstruments"),
                  secondaryPhenotype=c("Ramyg","Lamyg","Rhippo", "Lhippo", "Rpal", "Lpal", "Rput", "Lput"),
                  covariates=c("DEP_Z","Age_demeaned","Sex", "ICV_div1000_demeaned"),kinship=kinmat)

# test cortical surface area ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("Total_Zscore_bothinstruments"),
                  secondaryPhenotype=c("R_superiorfrontal_surfavg","L_superiorfrontal_surfavg","R_caudalmiddlefrontal_surfavg", "L_caudalmiddlefrontal_surfavg", "R_rostralmiddlefrontal_surfavg", "L_rostralmiddlefrontal_surfavg", "R_lateralorbitofrontal_surfavg", "L_lateralorbitofrontal_surfavg",
                                       "R_medialorbitofrontal_surfavg", "L_medialorbitofrontal_surfavg", "R_precentral_surfavg", "L_precentral_surfavg",
                                       "R_caudalanteriorcingulate_surfavg", "L_caudalanteriorcingulate_surfavg", "R_rostralanteriorcingulate_surfavg", "L_rostralanteriorcingulate_surfavg",
                                       "R_rostralanteriorcingulate_surfavg", "L_rostralanteriorcingulate_surfavg", "R_insula_surfavg", "L_insula_surfavg",
                                       "R_superiorparietal_surfavg", "L_superiorparietal_surfavg","R_inferiorparietal_surfavg", "L_inferiorparietal_surfavg",
                                       "R_precuneus_surfavg", "L_precuneus_surfavg", "R_supramarginal_surfavg", "L_supramarginal_surfavg",
                                       "R_postcentral_surfavg", "L_postcentral_surfavg", "R_temporalpole_surfavg", "L_temporalpole_surfavg", 
                                       "R_inferiortemporal_surfavg", "L_inferiortemporal_surfavg", "R_superiortemporal_surfavg", "L_superiortemporal_surfavg",
                                       "R_fusiform_surfavg", "L_fusiform_surfavg", "R_transversetemporal_surfavg", "L_transversetemporal_surfavg"),
                  covariates=c("DEP_Z","Age_demeaned","Sex", "TotalGlobalSurfArea_div1000_demeaned"),kinship=kinmat)

# test cortical thickness ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("Total_Zscore_bothinstruments"),
                  secondaryPhenotype=c("R_superiorfrontal_thickavg","L_superiorfrontal_thickavg","R_caudalmiddlefrontal_thickavg", "L_caudalmiddlefrontal_thickavg", "R_rostralmiddlefrontal_thickavg", "L_rostralmiddlefrontal_thickavg", "R_lateralorbitofrontal_thickavg", "L_lateralorbitofrontal_thickavg",
                                       "R_medialorbitofrontal_thickavg", "L_medialorbitofrontal_thickavg", "R_precentral_thickavg", "L_precentral_thickavg",
                                       "R_caudalanteriorcingulate_thickavg", "L_caudalanteriorcingulate_thickavg", "R_rostralanteriorcingulate_thickavg", "L_rostralanteriorcingulate_thickavg",
                                       "R_rostralanteriorcingulate_thickavg", "L_rostralanteriorcingulate_thickavg", "R_insula_thickavg", "L_insula_thickavg",
                                       "R_superiorparietal_thickavg", "L_superiorparietal_thickavg","R_inferiorparietal_thickavg", "L_inferiorparietal_thickavg",
                                       "R_precuneus_thickavg", "L_precuneus_thickavg", "R_supramarginal_thickavg", "L_supramarginal_thickavg",
                                       "R_postcentral_thickavg", "L_postcentral_thickavg", "R_temporalpole_thickavg", "L_temporalpole_thickavg", 
                                       "R_inferiortemporal_thickavg", "L_inferiortemporal_thickavg", "R_superiortemporal_thickavg", "L_superiortemporal_thickavg",
                                       "R_fusiform_thickavg", "L_fusiform_thickavg", "R_transversetemporal_thickavg", "L_transversetemporal_thickavg"),
                  covariates=c("DEP_Z","Age_demeaned","Sex", "GlobalMeanThickness_demeaned"),kinship=kinmat)

#### Everything with FNE
# test basic imaging phenotypes 
automatizedLmekin(dataset,primaryPhenotype=c("FNE"),secondaryPhenotype=c("TotalGlobalSurfArea_div1000_demeaned", "ICV_div1000_demeaned", "GlobalMeanThickness_demeaned"),covariates=c("Age_demeaned","Sex"),kinship=kinmat)

# test subcortical volumes ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("FNE"),
                  secondaryPhenotype=c("Ramyg","Lamyg","Rhippo", "Lhippo", "Rpal", "Lpal", "Rput", "Lput"),
                  covariates=c("DEP_Z","Age_demeaned","Sex", "ICV_div1000_demeaned"),kinship=kinmat)

# test cortical surface area ROIs 
automatizedLmekin(dataset,primaryPhenotype=c("FNE"),
                  secondaryPhenotype=c("R_superiorfrontal_surfavg","L_superiorfrontal_surfavg","R_caudalmiddlefrontal_surfavg", "L_caudalmiddlefrontal_surfavg", "R_rostralmiddlefrontal_surfavg", "L_rostralmiddlefrontal_surfavg", "R_lateralorbitofrontal_surfavg", "L_lateralorbitofrontal_surfavg",
                                       "R_medialorbitofrontal_surfavg", "L_medialorbitofrontal_surfavg", "R_precentral_surfavg", "L_precentral_surfavg",
                                       "R_caudalanteriorcingulate_surfavg", "L_caudalanteriorcingulate_surfavg", "R_rostralanteriorcingulate_surfavg", "L_rostralanteriorcingulate_surfavg",
                                       "R_rostralanteriorcingulate_surfavg", "L_rostralanteriorcingulate_surfavg", "R_insula_surfavg", "L_insula_surfavg",
                                       "R_superiorparietal_surfavg", "L_superiorparietal_surfavg","R_inferiorparietal_surfavg", "L_inferiorparietal_surfavg",
                                       "R_precuneus_surfavg", "L_precuneus_surfavg", "R_supramarginal_surfavg", "L_supramarginal_surfavg",
                                       "R_postcentral_surfavg", "L_postcentral_surfavg", "R_temporalpole_surfavg", "L_temporalpole_surfavg", 
                                       "R_inferiortemporal_surfavg", "L_inferiortemporal_surfavg", "R_superiortemporal_surfavg", "L_superiortemporal_surfavg",
                                       "R_fusiform_surfavg", "L_fusiform_surfavg", "R_transversetemporal_surfavg", "L_transversetemporal_surfavg"),
                  covariates=c("DEP_Z","Age_demeaned","Sex", "TotalGlobalSurfArea_div1000_demeaned"),kinship=kinmat)

# test cortical thickness ROIs
automatizedLmekin(dataset,primaryPhenotype=c("FNE"),
                  secondaryPhenotype=c("R_superiorfrontal_thickavg","L_superiorfrontal_thickavg","R_caudalmiddlefrontal_thickavg", "L_caudalmiddlefrontal_thickavg", "R_rostralmiddlefrontal_thickavg", "L_rostralmiddlefrontal_thickavg", "R_lateralorbitofrontal_thickavg", "L_lateralorbitofrontal_thickavg",
                                       "R_medialorbitofrontal_thickavg", "L_medialorbitofrontal_thickavg", "R_precentral_thickavg", "L_precentral_thickavg",
                                       "R_caudalanteriorcingulate_thickavg", "L_caudalanteriorcingulate_thickavg", "R_rostralanteriorcingulate_thickavg", "L_rostralanteriorcingulate_thickavg",
                                       "R_rostralanteriorcingulate_thickavg", "L_rostralanteriorcingulate_thickavg", "R_insula_thickavg", "L_insula_thickavg",
                                       "R_superiorparietal_thickavg", "L_superiorparietal_thickavg","R_inferiorparietal_thickavg", "L_inferiorparietal_thickavg",
                                       "R_precuneus_thickavg", "L_precuneus_thickavg", "R_supramarginal_thickavg", "L_supramarginal_thickavg",
                                       "R_postcentral_thickavg", "L_postcentral_thickavg", "R_temporalpole_thickavg", "L_temporalpole_thickavg", 
                                       "R_inferiortemporal_thickavg", "L_inferiortemporal_thickavg", "R_superiortemporal_thickavg", "L_superiortemporal_thickavg",
                                       "R_fusiform_thickavg", "L_fusiform_thickavg", "R_transversetemporal_thickavg", "L_transversetemporal_thickavg"),
                  covariates=c("DEP_Z","Age_demeaned","Sex", "GlobalMeanThickness_demeaned"),kinship=kinmat)

########################################################################################################################

