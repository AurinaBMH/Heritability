
#----------------------------------------------------------------------------------------------------------------------
#                                                             Twin simulate
# DVM Bishop, 11th March 2010, Based on script in OpenMXUserGuide, p 15
#----------------------------------------------------------------------------------------------------------------------
require(OpenMx)   # not needed for the simulation, but will be needed when we come to model specification
require(MASS)    # needed for multivariate random number generation
library(R.matlab)
set.seed(200)        # specified seed ensures same random number set generated on each run

parcellation = "HCPMMP1" #custom200" #"HCPMMP1"
tract = "iFOD2"
weights = "FA"

options(warn=1)
setwd("~/GoogleDrive/Genetics_connectome/Heritability/data/output")
data_covar = readMat("~/GoogleDrive/Genetics_connectome/Heritability/data/general/twinCovariatesDWI.mat")
data = readMat(sprintf("~/GoogleDrive/Genetics_connectome/Heritability/data/output/twinEdges_%s_%s_%sTEST5973.mat",parcellation, tract, weights))

# mark 0 values as NA because these will be the outliers in the data - we want to exclude them.
data$Output.DZ[data$Output.DZ == 0] <- NA
data$Output.MZ[data$Output.MZ == 0] <- NA

#data$Output.DZ <- data$Output.DZ/100
#data$Output.MZ <- data$Output.MZ/100
# all edges
numEdges = dim(data$Output.DZ)[3]
heritabilityA <- numeric(numEdges)
heritabilityC <- numeric(numEdges)
heritabilityE <- numeric(numEdges)

PvalsSat = data.frame(matrix(ncol = 3, nrow = numEdges))
PvalsACE = data.frame(matrix(ncol = 2, nrow = numEdges))
colnames(PvalsSat) <- c("ACE", "CE", "AE")
colnames(PvalsACE) <- c("CE", "AE")

mxOption(NULL,"Default optimizer","SLSQP")
#mxOption(NULL,"Default optimizer","NPSOL")
 
# change age from years to centuries
data_covar$DZ.age = data_covar$DZ.age/100
data_covar$MZ.age = data_covar$MZ.age/100

p_AE <- numeric(numEdges)
p_CE <- numeric(numEdges)
p_AC <- numeric(numEdges)

# fit models for each edge
for (edge in c(1:numEdges)){
  
myDataMZ<-data.frame(data$Output.MZ[,,edge], data_covar$MZ.age[,1], data_covar$MZ.age[,2], data_covar$MZ.sex[,1], data_covar$MZ.sex[,2])
myDataDZ<-data.frame(data$Output.DZ[,,edge], data_covar$DZ.age[,1], data_covar$DZ.age[,2], data_covar$DZ.sex[,1], data_covar$DZ.sex[,2])

myDataMZ_measure<-data$Output.MZ[,,edge]
myDataDZ_measure<-data$Output.DZ[,,edge]

colnames(myDataMZ) <- c('twin1', 'twin2','ageT1MZ', 'ageT2MZ', 'sexT1MZ', 'sexT2MZ')
colnames(myDataDZ) <- c('twin1', 'twin2','ageT1DZ', 'ageT2DZ', 'sexT1DZ', 'sexT2DZ')
selVars <- c('twin1','twin2')

# "complete" specifies use only cases with data in all columns
CovMZ = cov(data$Output.MZ[,,edge],use="complete")
CovDZ = cov(data$Output.DZ[,,edge],use="complete")

MeanMZ = colMeans(data$Output.MZ[,,edge],na.rm=TRUE)
MeanMZ = mean(MeanMZ)
MeanDZ = colMeans(data$Output.DZ[,,edge],na.rm=TRUE)
MeanDZ = mean(MeanDZ)


# do scatterplots for MZ and DZ
#split.screen(c(1,2))        # split display into two screens side by side
#screen(1)
#plot(myDataMZ_measure,main='MZ')    # main specifies overall plot title
#screen(2)
#plot(myDataDZ_measure,main='DZ')

#--------------------------------------------------------------------------------------------
#  Define saturated model
#--------------------------------------------------------------------------------------------
mylabels=c("twin1","twin2")
MZsat <- mxModel("MZsat",
                 
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(MeanMZ,MeanMZ), labels =c("b0_mz1","b0_mz2"), name="Intercepts" ),
                 mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                 mxMatrix( type="Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ"), name="MZDefVars"),
                 mxAlgebra( expression=Intercepts + beta %*% MZDefVars, name="expMeanMZ"),
                 mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, 0.5, name="CholMZ" ),
                 mxAlgebra( CholMZ %*% t(CholMZ), name="expCovMZ"),
                 mxData( myDataMZ, type="raw"), mxFitFunctionML(),
                 mxExpectationNormal( "expCovMZ", "expMeanMZ", mylabels))

DZsat <- mxModel("DZsat",
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(MeanDZ,MeanDZ), labels=c("b0_dz1","b0_dz2"), name="Intercepts" ),
                 mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                 mxMatrix( type="Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ"), name="DZDefVars"),
                 mxAlgebra( expression=Intercepts + beta %*% DZDefVars, name="expMeanDZ"),
                 mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, 0.5, name="CholDZ" ),
                 mxAlgebra( CholDZ %*% t(CholDZ),name="expCovDZ"),
                 mxData( myDataDZ, type="raw"), mxFitFunctionML(),
                 mxExpectationNormal( "expCovDZ", "expMeanDZ", mylabels))

mylabels=c("twin1","twin2")

# Model specification starts here
SatModel <- mxModel("twinSat", MZsat, DZsat, mxFitFunctionMultigroup(c('MZsat', 'DZsat')))
# adds together likelihoods for MZ and DZ groups
#---------------------------------------------------------------------------------------------------------------------------
#SatModelFit <- mxTryHard(SatModel) #The mxTryHard command evaluates the model.
#summary(SatModelFit)

#Equate intercepts across twin order:
sub1 <- omxSetParameters(model=SatModel,labels=c("b0_mz1","b0_mz2"),newlabels="b0_mz",name="Submodel1")
sub1 <- omxSetParameters(model=sub1,labels=c("b0_dz1","b0_dz2"),newlabels="b0_dz",name="Submodel1")
#sub1Fit <- mxRun(sub1)
#mxCompare(SatModelFit, sub1Fit) # compare models

#Equate intercepts across zygosity:
#sub2 <- omxSetParameters(model=sub1,labels=c("b0_mz","b0_dz"),newlabels="b0",name="Submodel2")
#sub2Fit <- mxTryHard(sub2)
#mxCompare(sub1Fit, sub2Fit) # compare models

#---------------------------------------------------------------------------------------------------------------------------
#this part doesn't work for the current parameterization of the covariance matrices
#Equate variance across twin order:
#sub3 <- omxSetParameters(model=sub2,labels=c("mzv1","mzv2"),newlabels="mzv",name="Submodel3")
#sub3 <- omxSetParameters(model=sub3,labels=c("dzv1","dzv2"),newlabels="dzv",name="Submodel3")

#Equate variance across zygosity:
#sub4 <- omxSetParameters(model=sub3,labels=c("mzv","dzv"),newlabels="v",name="Submodel4")

# -----------------------------------------------------------------------
#Fit ACE Model with RawData and Matrices Input
# -----------------------------------------------------------------------
twinACE <- mxModel("twinACE",
                   # Matrices X, Y, and Z to store a, c, and e path coefficients
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[2,2]/3)), label="a", name="X" ),
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[2,2]/3)), label="c", name="Y" ),
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[2,2]/3)), label="e", name="Z" ),
                   # Matrices A, C, and E compute variance components
                   mxAlgebra( expression=X %*% t(X), name="A" ),
                   mxAlgebra( expression=Y %*% t(Y), name="C" ),
                   mxAlgebra( expression=Z %*% t(Z), name="E" ),
                   
                   #mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values = MeanDZ, label="mean", name="expMean" ),
                   # Declare a matrix for the definition variable regression parameters, called beta
                   #mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, label=c("betaAge","betaSex"), name="beta"),
                   
                   # Algebra for expected variance/covariance matrix in MZ
                   mxAlgebra(
                     expression= rbind  (cbind(A+C+E , A+C),
                                         cbind(A+C   , A+C+E)),
                     name="expCovMZ"),
                   
                   # Algebra for expected variance/covariance matrix in DZ
                   # note use of 0.5, converted to 1*1 matrix
                   mxAlgebra(
                     expression= rbind  (cbind(A+C+E , 0.5%x%A+C),
                                         cbind(0.5%x%A+C , A+C+E)),
                     name="expCovDZ"),
                   
                   mxModel("MZ", mxData( observed=myDataMZ, type="raw" ),
                           # Algebra for making the means a function of the definition variables age and sex
                           mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(MeanMZ,MeanMZ), labels =c("b0_mz1","b0_mz2"), name="Intercepts" ),
                           #mxMatrix( type = "Full", nrow=1, ncol=1, free=T, c(MeanMZ), labels =c("b0"), name="Intercepts" ),
                           mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                           mxMatrix( type="Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ"), name="MZDefVars"),
                           mxAlgebra( expression=Intercepts + beta %*% MZDefVars, name="expMeanMZ"),
                           mxExpectationNormal( covariance="twinACE.expCovMZ", means="expMeanMZ", dimnames=selVars ),
                           mxFitFunctionML()),
                   
                   mxModel("DZ", mxData( observed=myDataDZ, type="raw" ),
                           mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(MeanDZ,MeanDZ), labels=c("b0_dz1","b0_dz2"), name="Intercepts" ),
                           #mxMatrix( type = "Full", nrow=1, ncol=1, free=T, c(MeanMZ), labels=c("b0"), name="Intercepts" ),
                           mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                           mxMatrix( type="Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ"), name="DZDefVars"),
                           mxAlgebra( expression=Intercepts + beta %*% DZDefVars, name="expMeanDZ"),
                           mxExpectationNormal( covariance="twinACE.expCovDZ", means="expMeanDZ", dimnames=selVars ),
                           mxFitFunctionML()),
                   
                   mxFitFunctionMultigroup( c("MZ.fitfunction",  "DZ.fitfunction"))
)
twinACEFit<-mxTryHard(twinACE)
#Run ACE model
# -----------------------------------------------------------------------
estCovMZ  <- mxEval(twinACE.expCovMZ, twinACEFit)      # expected covariance matrix for MZ's
estCovDZ  <- mxEval(twinACE.expCovDZ, twinACEFit)      # expected covariance matrix for DZ's
estVA     <- mxEval(a*a, twinACEFit)              # additive genetic variance, a^2
estVC     <- mxEval(c*c, twinACEFit)              # shared enviromnemtal variance, c^2
estVE     <- mxEval(e*e, twinACEFit)              # unique environmental variance, e^2
estVP     <- (estVA+estVC+estVE)                  # total variance
estPropVA <- estVA/estVP                          # standardized additive genetic variance
estPropVC <- estVC/estVP                          # standardized shared enviromnemtal variance
estPropVE <- estVE/estVP                          # standardized unique environmental variance
estACE    <- rbind(cbind(estVA,estVC,estVE),      # table of estimates
                   cbind(estPropVA,estPropVC,estPropVE))
LL_ACE    <- mxEval(objective, twinACEFit)        # likelihood of ADE model

heritabilityA[edge] <- estPropVA
heritabilityC[edge] <- estPropVC
heritabilityE[edge] <- estPropVE
#summary(twinACEFit)

# Generate AE Model - C=0

twinAE <- twinACE
twinAE   <- mxRename(twinAE, "twinAE")
twinAE   <- omxSetParameters(twinAE, labels="c", free=FALSE, values=0 )
#AEFit     <- mxTryHard(twinAE)
#AESumm   <- summary(AEFit)


#Generate CE Model - A=0
twinCE   <- twinACE
twinCE   <- mxRename(twinCE, "twinCE")
twinCE   <- omxSetParameters(twinCE, labels="a", free=FALSE, values=0 )
#CEFit     <- mxTryHard(twinCE)
#CESumm   <- summary(CEFit)


#Generate AC Model, E=0 # this model fails to run with an error ("All fit attempts resulted in errors - check starting values or model specification")
#twinAC   <- twinACE
#twinAC   <- mxRename(twinAC, "twinAC")
#twinAC   <- omxSetParameters(twinAC, labels="e", free=FALSE, values=0 )
#ACFit     <- mxTryHard(twinAC)
#ACSumm   <- summary(ACFit)
#mxCompare(twinACEFit, ACFit)

# model comparison
options('digits' = 5)
# compare saturated model to ACE model to see if ACE is significantlly worse that saturated
#compValuesSat = mxCompare(SatModelFit, twinACEFit)
# compare ACE to AC and CE (given EC failed) to see if a simpler model still fits the data not significantly worse than the full
#compValuesACE = mxCompare(twinACEFit, c(CEFit, AEFit))
# compare saturated model to all other models
#compValuesACEsat = mxCompare(SatModelFit, c(twinACEFit,CEFit, AEFit))
# save p and AIC values form the comparisons
#PvalsSat[edge,] = compValuesACEsat$p[2:4]
#PvalsACE[edge,] = compValuesACE$p[2:3]
#AICvals[edge,] = c(compValuesSat$AIC[2],compValuesACE$AIC[2:3])
}

heritabilityACE <- data.frame(heritabilityA,heritabilityC,heritabilityE)
write.csv(heritabilityACE,"heritabilityACE_TESTHCP5973.txt",row.names=FALSE)
