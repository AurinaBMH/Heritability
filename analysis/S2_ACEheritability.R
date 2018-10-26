#! /Library/Frameworks/R.framework/Resources/bin/Rscript
S2_ACEheritability <- function (parcellation, tract, weights, cvMeasure, conDens) {
  
require(OpenMx)   # not needed for the simulation, but will be needed when we come to model specification
require(MASS)    # needed for multivariate random number generation
library(R.matlab)
set.seed(200)    

setwd("~/GoogleDrive/Genetics_connectome/Heritability/data/output")
data_covar = readMat("~/GoogleDrive/Genetics_connectome/Heritability/data/general/twinCovariatesDWI.mat")
data = readMat(sprintf("~/GoogleDrive/Genetics_connectome/Heritability/data/output/twinEdges_%s_%s_%s_%s%d.mat",parcellation, tract, weights, cvMeasure, round(conDens*100)))

# remove second additional sibling from the data
data$Output.MZ <- data$Output.MZ[,-4,];
data$Output.DZ <- data$Output.DZ[,-4,]; 
 
data_covar$MZ.ID <- data_covar$MZ.ID[,-4]; 
data_covar$MZ.age <- data_covar$MZ.age[,-4]; 
data_covar$MZ.sex <- data_covar$MZ.sex[,-4]; 

data_covar$DZ.ID <- data_covar$DZ.ID[,-4]; 
data_covar$DZ.age <- data_covar$DZ.age[,-4]; 
data_covar$DZ.sex <- data_covar$DZ.sex[,-4]; 

# replace years with centuries
data_covar$DZ.age = data_covar$DZ.age/100
data_covar$MZ.age = data_covar$MZ.age/100

# replace missing covariates with a very different number setting them to "pseudo-missing" values
data_covar$MZ.ID[is.na(data_covar$MZ.ID)] <- -999
data_covar$MZ.age[is.na(data_covar$MZ.age)] <- -999
data_covar$MZ.sex[is.na(data_covar$MZ.sex)] <- -999

data_covar$DZ.ID[is.na(data_covar$DZ.ID)] <- -999
data_covar$DZ.age[is.na(data_covar$DZ.age)] <- -999
data_covar$DZ.sex[is.na(data_covar$DZ.sex)] <- -999

# mark 0 values as NA because these will be the outliers in the data - we want to exclude them.
data$Output.DZ[data$Output.DZ == 0] <- NA
data$Output.MZ[data$Output.MZ == 0] <- NA

# all edges
numEdges = dim(data$Output.DZ)[3]
heritabilityA <- numeric(numEdges)
heritabilityC <- numeric(numEdges)
heritabilityE <- numeric(numEdges)
heritabilityS <- numeric(numEdges)
#heritabilityT <- numeric(numEdges)

#PvalsSat = data.frame(matrix(ncol = 3, nrow = numEdges))
#PvalsACE = data.frame(matrix(ncol = 2, nrow = numEdges))
#colnames(PvalsSat) <- c("ACE", "CE", "AE")
#colnames(PvalsACE) <- c("CE", "AE")

#p_AE <- numeric(numEdges)
#p_CE <- numeric(numEdges)
#p_AC <- numeric(numEdges)

mxOption(NULL,"Default optimizer","SLSQP")
#mxOption(NULL,"Default optimizer","NPSOL")

for (edge in c(1:numEdges)){
#for (edge in c(1:5)){
# each connection (phenotype of interest) is denoted as edge
myDataMZ<-data.frame(data$Output.MZ[,,edge], data_covar$MZ.age[,1], data_covar$MZ.age[,2],data_covar$MZ.age[,3],data_covar$MZ.sex[,1],data_covar$MZ.sex[,2],data_covar$MZ.sex[,3])
myDataDZ<-data.frame(data$Output.DZ[,,edge], data_covar$DZ.age[,1], data_covar$DZ.age[,2],data_covar$DZ.age[,3],data_covar$DZ.sex[,1],data_covar$DZ.sex[,2],data_covar$DZ.sex[,3])

myDataMZ_measure<-data$Output.MZ[,,edge]
myDataDZ_measure<-data$Output.DZ[,,edge]

colnames(myDataMZ) <- c('twin1', 'twin2', 'sib','ageT1MZ', 'ageT2MZ', 'ageSIBMZ', 'sexT1MZ', 'sexT2MZ', 'sexSIBMZ')
colnames(myDataDZ) <- c('twin1', 'twin2', 'sib','ageT1DZ', 'ageT2DZ', 'ageSIBDZ', 'sexT1DZ', 'sexT2DZ', 'sexSIBDZ')
selVars <- c('twin1','twin2', 'sib')

# "complete" specifies use only cases with data in all columns
CovMZ = cov(data$Output.MZ[,,edge],use="complete")
CovDZ = cov(data$Output.DZ[,,edge],use="complete")

# mean across MZ twins
MeanMZ = colMeans(data$Output.MZ[,1:2,edge],na.rm=TRUE)
MeanMZ = mean(MeanMZ)
# mean across DZ twins
MeanDZ = colMeans(data$Output.DZ[,1:2,edge],na.rm=TRUE)
MeanDZ = mean(MeanDZ)
# mean across siblings
MeanSIBMZ = mean(data$Output.MZ[,3,edge],na.rm=TRUE)
MeanSIBDZ = mean(data$Output.DZ[,3,edge],na.rm=TRUE)

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
                   
                   # Algebra for expected variance/covariance matrix in MZ+sibling
                   mxAlgebra(
                     expression= rbind  (cbind(A+C+E , A+C, 0.5%x%A+C),
                                         cbind(A+C   , A+C+E, 0.5%x%A+C),
                                         cbind(0.5%x%A+C, 0.5%x%A+C, A+C+E)),
                     name="expCovMZ"),
                   
                   # Algebra for expected variance/covariance matrix in DZ+siblig
                   mxAlgebra(
                     expression= rbind  (cbind(A+C+E , 0.5%x%A+C, 0.5%x%A+C),
                                         cbind(0.5%x%A+C , A+C+E, 0.5%x%A+C),
                                         cbind(0.5%x%A+C , 0.5%x%A+C, A+C+E)),
                     name="expCovDZ"),
                   
                   mxModel("MZ", mxData( observed=myDataMZ, type="raw" ),
                           # Algebra for making the means a function of the definition variables age and sex
                           mxMatrix( type="Full", nrow=1, ncol=3, free=T, c(MeanMZ,MeanMZ,MeanSIBMZ), labels =c("b0_mz1","b0_mz2","b0_mzsib"), name="Intercepts"),
                           mxMatrix( type="Full", nrow=1, ncol=2, free=T, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                           mxMatrix( type="Full", nrow=2, ncol=3, free=F, labels=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ","data.ageSIBMZ","data.sexSIBMZ"), name="MZDefVars"),
                           mxAlgebra( expression=Intercepts + beta %*% MZDefVars, name="expMeanMZ"),
                           mxExpectationNormal( covariance="twinACE.expCovMZ", means="expMeanMZ", dimnames=selVars ),
                           mxFitFunctionML()),
                   
                   mxModel("DZ", mxData( observed=myDataDZ, type="raw" ),
                           mxMatrix( type="Full", nrow=1, ncol=3, free=T, c(MeanDZ,MeanDZ,MeanSIBDZ), labels=c("b0_dz1","b0_dz2","b0_dzsib"), name="Intercepts"),
                           mxMatrix( type="Full", nrow=1, ncol=2, free=T, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                           mxMatrix( type="Full", nrow=2, ncol=3, free=F, labels=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ","data.ageSIBDZ","data.sexSIBDZ"), name="DZDefVars"),
                           mxAlgebra( expression=Intercepts + beta %*% DZDefVars, name="expMeanDZ"),
                           mxExpectationNormal( covariance="twinACE.expCovDZ", means="expMeanDZ", dimnames=selVars ),
                           mxFitFunctionML()),
                   
                   mxFitFunctionMultigroup( c("MZ.fitfunction",  "DZ.fitfunction"))
)
twinACEFit<-mxTryHard(twinACE)

#Run ACE model
# -----------------------------------------------------------------------
#estCovMZ  <- mxEval(twinACE.expCovMZ, twinACEFit)      # expected covariance matrix for MZ's
#estCovDZ  <- mxEval(twinACE.expCovDZ, twinACEFit)      # expected covariance matrix for DZ's
estVA1     <- mxEval(a*a, twinACEFit)              # additive genetic variance, a^2
estVC1     <- mxEval(c*c, twinACEFit)              # shared enviromnemtal variance, c^2
estVE1     <- mxEval(e*e, twinACEFit)              # unique environmental variance, e^2
estVP1     <- (estVA1+estVC1+estVE1)                  # total variance

heritmodels <- matrix(nrow = 3, ncol = 4)

heritmodels[1,1] <- estVA1/estVP1                          # standardized additive genetic variance
heritmodels[1,2] <- estVC1/estVP1  
heritmodels[1,3] <- estVE1/estVP1
heritmodels[1,4] <- twinACEFit@output$status$code
# Generate AE Model - C=0
twinAE <- twinACE
twinAE   <- mxRename(twinAE, "twinAE")
twinAE   <- omxSetParameters(twinAE, labels="c", free=FALSE, values=0 )
twinAEFit    <- mxTryHard(twinAE)
AESumm   <- summary(twinAEFit)

estVA2    <- mxEval(a*a, twinAEFit)              # additive genetic variance, a^2
estVC2    <- mxEval(c*c, twinAEFit)              # shared enviromnemtal variance, c^2
estVE2    <- mxEval(e*e, twinAEFit)              # unique environmental variance, e^2
estVP2    <- (estVA2+estVC2+estVE2)                  # total variance

heritmodels[2,1] <- estVA2/estVP2                          # standardized additive genetic variance
heritmodels[2,2] <- estVC2/estVP2  
heritmodels[2,3] <- estVE2/estVP2
heritmodels[2,4] <- twinAEFit@output$status$code

#Generate CE Model - A=0
twinCE   <- twinACE
twinCE   <- mxRename(twinCE, "twinCE")
twinCE   <- omxSetParameters(twinCE, labels="a", free=FALSE, values=0 )
twinCEFit     <- mxTryHard(twinCE)

estVA3    <- mxEval(a*a, twinCEFit)              # additive genetic variance, a^2
estVC3    <- mxEval(c*c, twinCEFit)              # shared enviromnemtal variance, c^2
estVE3    <- mxEval(e*e, twinCEFit)              # unique environmental variance, e^2
estVP3    <- (estVA3+estVC3+estVE3)                  # total variance

heritmodels[3,1] <- estVA3/estVP3                          # standardized additive genetic variance
heritmodels[3,2] <- estVC3/estVP3  
heritmodels[3,3] <- estVE3/estVP3
heritmodels[3,4] <- twinCEFit@output$status$code

# model comparison
options('digits' = 5)
# compare ACE to AC and CE (given EC failed) to see if a simpler model still fits the data not significantly worse than the full
compValuesACE = mxCompare(twinACEFit, c(twinAEFit,twinCEFit))
# find model with the lowest AIC
INDmin = which.min(compValuesACE$AIC)
AICmin = compValuesACE$AIC[INDmin]

# get heritability estimates for the model with min AIC

heritabilityA[edge] <- heritmodels[INDmin,1]
heritabilityC[edge] <- heritmodels[INDmin,2]
heritabilityE[edge] <- heritmodels[INDmin,3]
heritabilityS[edge] <- heritmodels[INDmin,4]
}

heritabilityACE <- data.frame(heritabilityA,heritabilityC,heritabilityE,heritabilityS)
write.csv(heritabilityACE,sprintf("heritabilityACE_%s_%s_%s_%s%d.txt",parcellation, tract, weights, cvMeasure, round(conDens*100)),row.names=FALSE)
}