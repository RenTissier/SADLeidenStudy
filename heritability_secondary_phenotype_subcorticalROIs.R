

library(coxme)
library(mvtnorm)
library(lme4)
library(numDeriv)
library(car)
library(msm)
library(MASS)
library(kinship2)


Kinship.computation<-function(dataset,index){
  d2<-dataset[,index]
  ufam<-unique(as.vector(d2[,1]))
  nfam<-length(ufam)
  names(d2)<-c("fid", "iid", "pid", "mid", "sex")
  temp<-c(0,0,0,0,0)
  created<-rep(0,nfam)
  for (i in 1:nfam)
  { 
    pedi <- d2[d2[,1]==ufam[i],1:5]
    names(pedi)<-c("fid", "iid", "pid", "mid", "sex")
    fathers<-unique(pedi[,3])
    mothers<-unique(pedi[,4])
    father<-0
    mother<-0
    for (j in fathers)
      if (match(j, pedi[,2], nomatch = 0) == 0 && j!="0"){
        temp<-rbind(temp,data.frame(fid=ufam[i],iid=j,pid=0,mid=0,sex=1,stringsAsFactors=FALSE))
        father<-father+1 }
    for (k in mothers)
      if (match(k, pedi[,2], nomatch = 0) == 0 && k!="0"){
        temp<-rbind(temp,data.frame(fid=ufam[i],iid=k,pid=0,mid=0,sex=2,stringsAsFactors=FALSE))
        mother<-mother+1 }
    temp<-rbind(temp,pedi)
    created[i]<-mother+father
  }
  temp<-temp[-c(1),]
  I <- dim(temp)[1]
  a<-pedigree(id=temp$iid,dadid=temp$pid,momid=temp$mid,sex=temp$sex,famid=temp$fid)
  varFamList<-vector(mode="list",length=length(family))
  for (i in 1:nfam)
  {
    ped2 <- a[i]
    kinmat<-kinship(ped2)
    kinmat<-2*(kinmat[-c(1:created[i]),-c(1:created[i])])
    varFamList[[i]]<-kinmat
  }
  return(varFamList)
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


#########################################################################################################


covPhenotypes<-function(kinship,varianceA,varianceB,coeff,rho){
  size<-length(kinship[,1])	
  matrix<-rho*varianceA*varianceB*kinship
  matrix<-matrix+coeff^2*diag(size)
  matrix
}


Variance<-function(kinship,genetic,error,u){
  size<-length(kinship[,1])	
  matrix<-kinship*genetic^2+(error^2+u^2)*diag(size)
  matrix
}



comparison<-function(x,vec){
  a<-all(x==b)
  a
}



#######################################################################################################################
####################################################################################################################################################
CatalanoLikelihood<-function(dataset,VarX,VarY,covYX,sizeFamily){
  lower<-NULL
  upper<-NULL
  for(j in 1:sizeFamily){
    if(dataset$group[j]==1){
      lower<-c(lower,0)
      upper<-c(upper,+Inf)
      
    }
    if(dataset$group[j]==0){
      lower<-c(lower,-Inf)
      upper<-c(upper,0)
    }
  }	
  meanYgivenX<-as.numeric(dataset$meanY+covYX%*%solve(VarX)%*%(dataset$secondary-dataset$meanX)) #compute the mean of Y given the secondary phenotype
  VarYgivenX<-VarY-covYX%*%solve(VarX)%*%t(covYX) #compute the varicance-covariance matrix of Y|X
  prob<-dmvnorm(dataset$secondary,mean=dataset$meanX,sigma=VarX)*pmvnorm(lower,upper,mean=meanYgivenX,sigma=VarYgivenX,algorithm=GenzBretz(maxpts = 300000000,abseps =0.000000001))#Obtention of the loglikelihood for the family by summing log(P(X|G))+log(P(Y|G,X)))
  prob
}


denominatorBis<-function(dataset,VarY,sizeFamily){
  lower<-NULL
  upper<-NULL
  for (i in 1:sizeFamily){
    if (dataset$group[i]==1){
      lower<-c(lower,0)
      upper<-c(upper,+Inf)
    }
    else{
      lower<-c(lower,-Inf)
      upper<-c(upper,0)
    }
  }
  ProbY<-pmvnorm(lower,upper,mean=dataset$meanY,sigma=VarY,algorithm=GenzBretz(maxpts = 300000000,abseps =0.000000001)) #computing the marginal probability of Y given the genotypes
  ProbY
}
?pmvnorm

likelihoodFamily<-function(x,kinships,covariates,SecondaryPhenotype,PrimaryPhenotype,parameters){
  sizeFamily<-dim(kinships[[x]])[1]
  numberFixedEffects<-(length(covariates[[x]][1,])+1)
  parametersX<-parameters[1:(numberFixedEffects+2)]
  parametersX[(numberFixedEffects+1):(numberFixedEffects+2)]<-exp(parametersX[(numberFixedEffects+1):(numberFixedEffects+2)])
  parametersY<-parameters[(numberFixedEffects+3):(numberFixedEffects+4)]
  parametersY[2]<-exp(parametersY[2])
  parameterRho<-parameters[numberFixedEffects+5]
  parameterRho<-1/(1+exp(parameterRho))
  parameterXY<-parameters[numberFixedEffects+6]
  parameterXY<-exp(parameterXY)
  meanX<-as.matrix(cbind(rep(1,sizeFamily),covariates[[x]]))%*%parametersX[1:numberFixedEffects]
  meanY<-rep(parametersY[1],sizeFamily)
  VarX<-Variance(kinships[[x]],genetic=parametersX[(numberFixedEffects+1)],error=parametersX[(numberFixedEffects+2)],u=parameterXY)
  VarY<-Variance(kinships[[x]],genetic=parametersY[2],error=1,u=parameterXY)
  covYX<-covPhenotypes(kinships[[x]],varianceA=parametersX[(numberFixedEffects+1)],varianceB=parametersY[2],coeff=parameterXY,rho=parameterRho)
  b<-as.data.frame(cbind(meanX,meanY,group=PrimaryPhenotype[[x]],secondary=SecondaryPhenotype[[x]]))
  names(b)<-c("meanX","meanY","group","secondary")
  jointDensity<-CatalanoLikelihood(b,VarX,VarY,covYX,sizeFamily) #Computation of the joint probability
  numerator<-jointDensity
  denom<-denominatorBis(b,VarY,sizeFamily)
  resultsFamily<-log(numerator/denom)
  resultsFamily
}

likelihoodFamily0<-function(x,kinships,covariates,SecondaryPhenotype,PrimaryPhenotype,parameters){
  sizeFamily<-dim(kinships[[x]])[1]
  numberFixedEffects<-(length(covariates[[x]][1,])+1)
  parametersX<-parameters[1:(numberFixedEffects+1)]
  parametersX[(numberFixedEffects+1)]<-exp(parametersX[(numberFixedEffects+1)])
  parametersY<-parameters[(numberFixedEffects+2):(numberFixedEffects+3)]
  parametersY[2]<-exp(parametersY[2])
  parameterRho<-parameters[numberFixedEffects+4]
  parameterRho<-1/(1+exp(parameterRho))
  parameterXY<-parameters[numberFixedEffects+5]
  parameterXY<-exp(parameterXY)
  meanX<-as.matrix(cbind(rep(1,sizeFamily),covariates[[x]]))%*%parametersX[1:numberFixedEffects]
  meanY<-rep(parametersY[1],sizeFamily)
  VarX<-Variance(kinships[[x]],genetic=0,error=parametersX[(numberFixedEffects+1)],u=parameterXY)
  VarY<-Variance(kinships[[x]],genetic=parametersY[2],error=1,u=parameterXY)
  covYX<-covPhenotypes(kinships[[x]],varianceA=0,varianceB=parametersY[2],coeff=parameterXY,rho=parameterRho)
  b<-as.data.frame(cbind(meanX,meanY,group=PrimaryPhenotype[[x]],secondary=SecondaryPhenotype[[x]]))
  names(b)<-c("meanX","meanY","group","secondary")
  jointDensity<-CatalanoLikelihood(b,VarX,VarY,covYX,sizeFamily) #Computation of the joint probability
  numerator<-jointDensity
  denom<-denominatorBis(b,VarY,sizeFamily)
  resultsFamily<-log(numerator/denom)
  resultsFamily
}


LikelihoodComputation<-function(parameters,kinships,covariates,secondaryPhenotype,primaryPhenotype,seed){
  set.seed(seed) 
  a<-lapply(names(kinships),likelihoodFamily,kinships=kinships,PrimaryPhenotype=primaryPhenotype,SecondaryPhenotype=secondaryPhenotype,covariates=covariates,parameters=parameters)
  resa<-Reduce("+",a)
  resa
}

LikelihoodComputation0<-function(parameters,kinships,covariates,secondaryPhenotype,primaryPhenotype,seed){
  set.seed(seed) 
  a<-lapply(names(kinships),likelihoodFamily0,kinships=kinships,PrimaryPhenotype=primaryPhenotype,SecondaryPhenotype=secondaryPhenotype,covariates=covariates,parameters=parameters)
  resa<-Reduce("+",a)
  resa
}

LikelihoodOptimization<-function(dataset,primaryPhenotype,secondaryPhenotype,covariates,index,seed){
  set.seed(seed) 
  results<-vector("list", length = (length(secondaryPhenotype)+2))
  seeds<-vector("list", length = (length(secondaryPhenotype)+2))
  resheritability<-NULL
  negative<-FALSE

  for (j in 1:length(secondaryPhenotype)){
    variables<-as.vector(as.numeric(unique(c(index,which(names(dataset)%in%primaryPhenotype),which(names(dataset)%in%secondaryPhenotype[j]),which(names(dataset)%in%covariates)))))
    datasetPhenotype<-dataset[,variables]
    datasetPhenotype<-datasetPhenotype[complete.cases(datasetPhenotype),]
    datasetPhenotype[,secondaryPhenotype[j]]<-datasetPhenotype[,secondaryPhenotype[j]]/sd(datasetPhenotype[,secondaryPhenotype[j]])
    if (cor(datasetPhenotype[,6],datasetPhenotype[,7])<0){
      negative<-TRUE
      datasetPhenotype[,7]<--datasetPhenotype[,7]
    }

    varFamList<-Kinship.computation(datasetPhenotype,c(1:5))
    names(datasetPhenotype)[1]<-"fid"
    family<-unique(datasetPhenotype[,1])
    PrimaryPhenotype<-vector(mode="list",length=length(family))
    SecondaryPhenotype<-vector(mode="list",length=length(family))
    Covariates<-vector(mode="list",length=length(family))
    for (i in 1:length(family)){
      dfamily<-datasetPhenotype[datasetPhenotype[,1]==family[i],]
      PrimaryPhenotype[[i]]<-dfamily[,which(names(dfamily)%in%primaryPhenotype)]
      SecondaryPhenotype[[i]]<-dfamily[,which(names(dfamily)%in%secondaryPhenotype[j])]
      Covariates[[i]]<-as.data.frame(dfamily[,which(names(dfamily)%in%covariates)])
    }	
    names(Covariates)<-c(1:length(family))
    names(varFamList)<-c(1:length(family))
    names(PrimaryPhenotype)<-c(1:length(family))
    names(SecondaryPhenotype)<-c(1:length(family))
    
    Formula1 <- formula(paste(primaryPhenotype," ~  (1|fid)"))
    Formula2 <- formula(paste(secondaryPhenotype[j]," ~ ", 
                              paste(covariates, collapse=" + ")," + (1|fid)"))
    e<-glmer(Formula1,data=datasetPhenotype,family=binomial(link="probit"))
    coefs2 <- data.frame(coef(summary(e)))
    primaryFixed<-as.numeric(fixef(e))
    vcPrimary<-VarCorr(e)[[1]]
    f<-lmer(Formula2,data=datasetPhenotype)
    coefs <- data.frame(coef(summary(f)))
    secondaryFixed<-as.numeric(fixef(f))
    secondaryFixed<-secondaryFixed
    vcSecondary<-as.numeric(VarCorr(f)[[1]])
    residuals<-attr(VarCorr(f),"sc")
    standardErrorHeritability2<-NA
    pvalHeritability<-1
    heritability<-0
    seedPhenotype<-c(seed)
    print(c(paste("Analyzing phenotype ",secondaryPhenotype,sep=" ")))
    while(pvalHeritability==1){
      res0<-try(optim(c(secondaryFixed[1],rep(0,length(covariates)),log(residuals),primaryFixed[1],log(1),0,log(sqrt(2))),LikelihoodComputation0,kinships=varFamList,covariates=Covariates,secondaryPhenotype=SecondaryPhenotype,primaryPhenotype=PrimaryPhenotype,seed=seed,method="BFGS",control=list(abstol=1e-16,trace=6,fnscale=-1,parscale=c(rep(1,(length(covariates)+1)),0.1,1,0.1,0.1,0.1),maxit=200),hessian=T))
      res<-try(optim(c(secondaryFixed[1],rep(0,length(covariates)),log(1),log(residuals),primaryFixed[1],log(1),0,log(sqrt(2))),LikelihoodComputation,kinships=varFamList,covariates=Covariates,secondaryPhenotype=SecondaryPhenotype,primaryPhenotype=PrimaryPhenotype,seed=seed,method="BFGS",control=list(abstol=1e-16,trace=6,fnscale=-1,parscale=c(rep(1,(length(covariates)+1)),0.1,0.1,1,0.1,0.1,0.1),maxit=200),hessian=T))
      pvalHeritability<-1-pchisq(2*(res$value-res0$value),1)
      suppressWarnings(try(StdError<-sqrt(diag(ginv(-res$hessian)))))
      heritability<-exp(res$par[(length(covariates)+2)])^2/(exp(res$par[(length(covariates)+2)])^2+exp(res$par[(length(covariates)+3)])^2+exp(res$par[(length(covariates)+7)])^2)
      heritabilityPrimary<-exp(res$par[(length(covariates)+5)])^2/(exp(res$par[(length(covariates)+5)])^2+1+exp(res$par[(length(covariates)+7)])^2)
      try(covarianceMatrix<-ginv(-res$hessian))
      CovHeritability<-as.matrix(rbind(c(covarianceMatrix[(length(covariates)+2),(length(covariates)+2)],covarianceMatrix[(length(covariates)+2),(length(covariates)+3)],covarianceMatrix[(length(covariates)+2),(length(covariates)+7)]),c(covarianceMatrix[(length(covariates)+3),(length(covariates)+2)],covarianceMatrix[(length(covariates)+3),(length(covariates)+3)],covarianceMatrix[(length(covariates)+3),(length(covariates)+7)]),c(covarianceMatrix[(length(covariates)+7),(length(covariates)+2)],covarianceMatrix[(length(covariates)+7),(length(covariates)+3)],covarianceMatrix[(length(covariates)+7),(length(covariates)+7)])))
      standardErrorHeritability<-suppressWarnings(deltamethod(~(exp(x1)^2/(exp(x1)^2+exp(x2)^2+exp(x3)^2)),c(res$par[(length(covariates)+2)],res$par[(length(covariates)+3)],exp(res$par[(length(covariates)+7)])),CovHeritability))
      standardErrorHeritability2<-StdError[(length(covariates)+2)]
      #pvalHeritability<-1-pchisq(heritability^2/standardErrorHeritability^2,1)
      #pvalHeritability2<-1-pchisq(res$par[(length(covariates)+2)]^2/StdError[(length(covariates)+2)]^2,1)
      #if (!is.finite(pvalHeritability)){
      #  pvalHeritability<-0
      #}
      #if (!is.finite(standardErrorHeritability)){
      #  standardErrorHeritability<-0
      #}
      seed<-sample(1:100000,1)
      seedPhenotype<-c(seedPhenotype,seed)
    }
    numberFixedEffects<-(length(Covariates[[1]][1,])+1)
    if(negative){
      res$par[1:(numberFixedEffects)]<--res$par[1:(numberFixedEffects)]
      par<-c(res$par[1:(numberFixedEffects)],exp(res$par[(numberFixedEffects+1):(numberFixedEffects+2)]),res$par[(numberFixedEffects+3)],exp(res$par[(numberFixedEffects+4)]),1/(1+exp(res$par[(numberFixedEffects+5)])),exp(res$par[(numberFixedEffects+6)]))
    }
    if(!negative){
      par<-c(res$par[1:(numberFixedEffects)],exp(res$par[(numberFixedEffects+1):(numberFixedEffects+2)]),res$par[(numberFixedEffects+3)],exp(res$par[(numberFixedEffects+4)]),1/(1+exp(res$par[(numberFixedEffects+5)])),exp(res$par[(numberFixedEffects+6)]))
    }
    resOnePheno<-c(heritability,pvalHeritability)
    resheritability<-rbind(resheritability,resOnePheno)
    resCovariates<-NULL
    for(k in 1:(length(covariates))){
      standard<-StdError[(k+1)]
      value<-res$par[(k+1)] 
      pval<-1-pchisq(value^2/standard^2,1)
      covaRes<-c(value,standard,pval)
      resCovariates<-rbind(resCovariates,covaRes)
    }
    row.names(resCovariates)<-covariates
    colnames(resCovariates)<-c("estimate","standard error","p-value")
    results[[(j+1)]]<-resCovariates
    seeds[[j]]<-seedPhenotype
    write.table(t(c(secondaryPhenotype[j],resOnePheno,seedPhenotype[(length(seedPhenotype)-1)])),"heritabilityResults.txt",sep="\t",append=T,col.names=F)
    write.table(resCovariates,paste(secondaryPhenotype[j],"covariatesResults.txt",sep=""),sep="\t")
  }
  row.names(resheritability)<-secondaryPhenotype
  colnames(resheritability)<-c("heritability","p-value")
  print(resheritability)
  seedsPhenotype<-NULL
  for (k in 1:length(seeds)){
    seedsPhenotype<-c(seedsPhenotype,seeds[[k]][length(seeds[[k]])])
  }
  results[[1]]<-resheritability
  results[[(length(secondaryPhenotype)+2)]]<-seedsPhenotype
}

##########################################################Example###########################################################################################
#########importation of the data

dataset<-read.table("BasHoogendam_2018_SMRI_GM_dataset.txt",sep="\t",header=T) #importation of the data
dataset$Age_demeaned<-dataset$Age_demeaned/sd(dataset$Age_demeaned,na.rm=T)
dataset$ICV_div1000_demeaned<-dataset$ICV_div1000_demeaned/sd(dataset$ICV_div1000_demeaned,na.rm=T)
dataset$GlobalMeanThickness_demeaned<-dataset$GlobalMeanThickness_demeaned/sd(dataset$GlobalMeanThickness_demeaned,na.rm=T)
dataset$TotalGlobalSurfArea_div1000_demeaned<-dataset$TotalGlobalSurfArea_div1000_demeaned/sd(dataset$TotalGlobalSurfArea_div1000_demeaned,na.rm=T)
Phenotypes<-c("Lput", "Rput", "Lamyg", "Ramyg", "Lhippo", "Rhippo", "Lpal", "Rpal") # input all secondary phenotypes here

########function to use LikelihoodOptimization##################################################
######parameters:
######dataset= dataset to use
######primaryPhenotype = name of the case-control variable 
######secondaryPhenotype = name of the variable you want heritability estimates
######covariates = names of variable you want to adjust for
######index = vector of index containing IN THE GOOD ORDER the column number of Id family, Id individual, Id father, Id mother and Sex (same as for testing association)
# dataset$ICV_div1000<-dataset$ICV/1000

########Example with xx as secondary phenotype and Age and Sex and ICV as covariates, SAD_d as primary phenotype
results_subcortical<-LikelihoodOptimization(dataset,primaryPhenotype=c("SAD_d"),secondaryPhenotype=Phenotypes[m],covariates=c("Age_demeaned", "Sex", "ICV_div1000_demeaned"),index=c(2,3,5,4,7),seed=32965)


