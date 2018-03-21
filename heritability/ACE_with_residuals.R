install.packages("OpenMx")
install.packages("XLConnect")
require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
library(XLConnect)
library(OpenMx)
mxOption(NULL,"Default optimizer","SLSQP")
library(R.matlab)
options(warn=1)
setwd("~/GoogleDrive/Genetics_connectome/Heritability/data/general")
data_covar = readMat("twinCovariatesDWI.mat")
data = readMat("~/GoogleDrive/Genetics_connectome/Heritability/data/output/twinEdges_HCPMMP1_iFOD2_FA.mat")
# all edges
numEdges = dim(data$Output.DZ)[3]
heritabilityA <- numeric(numEdges)
heritabilityC <- numeric(numEdges)
heritabilityE <- numeric(numEdges)

for (edge in c(1:20)){

MeasureMZ = 'Edge weight MZ'
MeasureDZ = 'Edge weight DZ'

mzData_measure<-data$Output.MZ[,,edge]
dzData_measure<-data$Output.DZ[,,edge]

colnames(mzData_measure) <- c('twin1', 'twin2')
colnames(dzData_measure) <- c('twin1', 'twin2')

# Generate Descriptive Statistics means and covariances
MeanMZ = colMeans(mzData_measure,na.rm=TRUE)
MeanMZ = mean(MeanMZ)
MeanDZ = colMeans(dzData_measure,na.rm=TRUE)
MeanDZ = mean(MeanDZ)
CovMZ = cov(mzData_measure,use="complete")
CovDZ = cov(dzData_measure,use="complete")

# plot the data

pictureName = sprintf("~/GoogleDrive/Genetics_connectome/Heritability/data/output/%sedge.png",edge)
png (pictureName)
split.screen(c(1,2))
screen(1)
plot(mzData_measure,main=MeasureMZ,xlim=c(min(mzData_measure, dzData_measure), max(mzData_measure, dzData_measure)), ylim=c(min(mzData_measure, dzData_measure), max(mzData_measure, dzData_measure)))
screen(2)
plot(dzData_measure,main=MeasureDZ,xlim=c(min(mzData_measure, dzData_measure), max(mzData_measure, dzData_measure)), ylim=c(min(mzData_measure, dzData_measure), max(mzData_measure, dzData_measure)))
dev.off() # to complete the writing process and return output to your monitor
# -----------------------------------------
mzData<-data.frame(data$Output.MZ[,,edge], data_covar$MZ.age[,1], data_covar$MZ.age[,2], data_covar$MZ.sex[,1], data_covar$MZ.sex[,2])
dzData<-data.frame(data$Output.DZ[,,edge], data_covar$DZ.age[,1], data_covar$DZ.age[,2], data_covar$DZ.sex[,1], data_covar$DZ.sex[,2])

colnames(mzData) <- c('twin1', 'twin2', 'ageT1MZ', 'ageT2MZ', 'sexT1MZ', 'sexT2MZ')
colnames(dzData) <- c('twin1', 'twin2', 'ageT1DZ', 'ageT2DZ', 'sexT1DZ', 'sexT2DZ')
selVars <- c('twin1','twin2')

#Fit ACE Model with RawData and Matrices Input
# -----------------------------------------------------------------------
twinACE <- mxModel("twinACE",
                   # Matrices X, Y, and Z to store a, c, and e path coefficients
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[1,2]/3)), label="a", name="X" ),
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[1,2]/3)), label="c", name="Y" ),
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[1,2]/3)), label="e", name="Z" ),
                   # Matrices A, C, and E compute variance components
                   mxAlgebra( expression=X %*% t(X), name="A" ),
                   mxAlgebra( expression=Y %*% t(Y), name="C" ),
                   mxAlgebra( expression=Z %*% t(Z), name="E" ),

                   mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= MeanDZ, label="mean", name="expMean" ),
                   # Declare a matrix for the definition variable regression parameters, called beta
                   mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, label=c("betaAge","betaSex"), name="beta"),

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

                   mxModel("MZ", mxData( observed=mzData, type="raw" ),
                           # Algebra for making the means a function of the definition variables age and sex
                           mxMatrix( type="Full", nrow=2, ncol=2, free=F, label=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ"), name="MZDefVars"),
                           mxAlgebra( expression=twinACE.expMean + twinACE.beta %*% MZDefVars, name="expMeanMZ"),
                           mxExpectationNormal(covariance="twinACE.expCovMZ", means="expMeanMZ", dimnames=selVars ),
                           mxFitFunctionML()),

                   mxModel("DZ", mxData( observed=dzData, type="raw" ),
                           mxMatrix( type="Full", nrow=2, ncol=2, free=F, label=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ"), name="DZDefVars"),
                           mxAlgebra( expression=twinACE.expMean + twinACE.beta %*% DZDefVars, name="expMeanDZ"),
                           mxExpectationNormal( covariance="twinACE.expCovDZ", means="expMeanDZ", dimnames=selVars ),
                           mxFitFunctionML()),

                   mxFitFunctionMultigroup( c("MZ.fitfunction",  "DZ.fitfunction"))

)

#Run ACE model
# -----------------------------------------------------------------------
twinACEFit <- mxRun(twinACE)
# running the same model twice - increases heritability for the 2nd edge
twinACEFit <- mxRun(twinACEFit)


summary(twinACEFit)

# Generate ACE Model Output

estMean   <- mxEval(mean, twinACEFit)             # expected mean
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

mzresids  <- matrix(NA,nrow=nrow(mzData),ncol=2)
for(i in 1:nrow(mzData)){
  mzresids[i,] <- as.matrix(mzData[i,c("twin1","twin2")]) - as.matrix(mxEval(expression=twinACE.expMean + twinACE.beta %*% MZ.MZDefVars, model = twinACEFit, compute=T,defvar.row=i))
}
dzresids  <- matrix(NA,nrow=nrow(dzData),ncol=2)
for(i in 1:nrow(dzData)){
  dzresids[i,] <- as.matrix(dzData[i,c("twin1","twin2")]) - as.matrix(mxEval(expression=twinACE.expMean + twinACE.beta %*% DZ.DZDefVars, model = twinACEFit, compute=T,defvar.row=i))
}

pictureNameResiduals = sprintf("~/GoogleDrive/Genetics_connectome/Heritability/data/output/%sedge_residuals.png",edge)
png (pictureNameResiduals)

split.screen(c(1,2))
screen(1)
plot(mzresids,main="mzresids", xlim=c(min(mzresids, dzresids), max(mzresids, dzresids)), ylim=c(min(mzresids, dzresids), max(mzresids, dzresids)))
screen(2)
plot(dzresids,main="dzresids", xlim=c(min(mzresids, dzresids), max(mzresids, dzresids)), ylim=c(min(mzresids, dzresids), max(mzresids, dzresids)))
dev.off() # to complete the writing process and return output to your monitor

}
