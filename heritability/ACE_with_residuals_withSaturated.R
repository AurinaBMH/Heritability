install.packages("OpenMx")
install.packages("XLConnect")
require(OpenMx)
require(psych)

#Prepare Data
# -----------------------------------------------------------------------
library(XLConnect)
library(OpenMx)
#mxOption(NULL,"Default optimizer","SLSQP")
library(R.matlab)

# choose options
parcellation = "HCPMMP1" #"HCPMMP1"
tract = "iFOD2"
weights = "standard"

options(warn=1)
setwd("~/GoogleDrive/Genetics_connectome/Heritability/data/output")
data_covar = readMat("~/GoogleDrive/Genetics_connectome/Heritability/data/general/twinCovariatesDWI.mat")
data = readMat(sprintf("~/GoogleDrive/Genetics_connectome/Heritability/data/output/twinEdges_%s_%s_%s.mat",parcellation, tract, weights))
# all edges
numEdges = dim(data$Output.DZ)[3]
heritabilityA <- numeric(numEdges)
heritabilityC <- numeric(numEdges)
heritabilityE <- numeric(numEdges)

Pvals = data.frame(matrix(ncol = 2, nrow = numEdges))
colnames(Pvals) <- c("CE", "AE")

p_AE <- numeric(numEdges)
p_CE <- numeric(numEdges)
p_AC <- numeric(numEdges)

for (edge in c(5)){

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

pictureName = sprintf("~/GoogleDrive/Genetics_connectome/Heritability/data/output/plots/%sedge_%s_%s_%s.png",edge, parcellation, tract, weights)
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

# Fit saturated model where relationships between twin pairs are ignored - models for MZ and DZ twins are identical.
#--------------------------------------------------------------------------------------------
# Model specification starts here
mylabels=c("twin1","twin2")
#fitFunction <- mxFitFunctionML(rowDiagnostics=TRUE)
mxFitFunctionMultigroup(c("MZsat","DZsat"))
#mgFitFun <- 
fitFunction <- mxFitFunctionML(rowDiagnostics=TRUE)

MZsat <- mxModel(
  "MZsat",
  mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(0,0), labels=c("b0_mz1","b0_mz2"), name="Intercepts" ),
  mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, label=c("betaAge","betaSex"), name="beta"),
  mxMatrix( type="Full", nrow=2, ncol=2, free=F, label=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ"), name="MZDefVars"),
  mxAlgebra( expression=Intercepts + beta %*% MZDefVars, name="expMeanMZ"),
  mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, .5, name="CholMZ" ),
  mxAlgebra( CholMZ %*% t(CholMZ), name="expCovMZ"),
  mxData( mzData, type="raw"), mxFitFunctionML(),
  mxExpectationNormal( "expCovMZ", "expMeanMZ", mylabels))

DZsat <- mxModel(
  "DZsat",
  mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(0,0), labels=c("b0_dz1","b0_dz2"), name="Intercepts" ),
  mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, label=c("betaAge","betaSex"), name="beta"),
  mxMatrix( type="Full", nrow=2, ncol=2, free=F, label=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ"), name="DZDefVars"),
  mxAlgebra( expression=Intercepts + beta %*% DZDefVars, name="expMeanDZ"),
  mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, .5, name="CholDZ" ),
  mxAlgebra( CholDZ %*% t(CholDZ), name="expCovDZ"),
  mxData( dzData, type="raw"), mxFitFunctionML(),
  mxExpectationNormal( "expCovDZ", "expMeanDZ", mylabels))

mytwinSatModel <- mxModel(model='twinSat', MZsat, DZsat, mxFitFunctionMultigroup(c('MZsat', 'DZsat')))

saturatedOut <- mxTryHard(mytwinSatModel)
summary(saturatedOut)

#Fit ACE Model with RawData and Matrices Input
# -----------------------------------------------------------------------
twinACE <- mxModel("twinACE",
                   # Matrices X, Y, and Z to store a, c, and e path coefficients
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[1,1]/3)), label="a", name="X" ),
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[1,1]/3)), label="c", name="Y" ),
                   mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[1,1]/3)), label="e", name="Z" ),
                   # Matrices A, C, and E compute variance components
                   mxAlgebra( expression=X %*% t(X), name="A" ),
                   mxAlgebra( expression=Y %*% t(Y), name="C" ),
                   mxAlgebra( expression=Z %*% t(Z), name="E" ),

                   mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=MeanDZ, label="mean", name="expMean" ),
                   # Declare a matrix for the definition variable regression parameters, called beta
                   mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=0, label=c("betaAge","betaSex"), name="beta"),

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
twinACEFit <- mxTryHard(twinACE)
# running the same model twice - increases heritability for the 2nd edge
twinACEFit <- mxTryHard(twinACEFit)

summary(twinACEFit)

# compare to saturated model
compValues = mxCompare(saturatedOut, twinACEFit)

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

msize = nrow(mzData_measure)*ncol(mzData_measure)
dsize = nrow(dzData_measure)*ncol(dzData_measure)
DF_ACE = msize+dsize-nrow(twinACEFit@output$standardErrors)

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

pictureNameResiduals = sprintf("~/GoogleDrive/Genetics_connectome/Heritability/data/output/plots/%sedge_residuals_%s_%s_%s.png",edge, parcellation, tract, weights)
png (pictureNameResiduals)

split.screen(c(1,2))
screen(1)
plot(mzresids,main="mzresids", xlim=c(min(mzresids, dzresids), max(mzresids, dzresids)), ylim=c(min(mzresids, dzresids), max(mzresids, dzresids)))
screen(2)
plot(dzresids,main="dzresids", xlim=c(min(mzresids, dzresids), max(mzresids, dzresids)), ylim=c(min(mzresids, dzresids), max(mzresids, dzresids)))
dev.off() # to complete the writing process and return output to your monitor

# Generate AE Model - C=0
twinAE   <- twinACE
twinAE   <- omxSetParameters( twinAE, labels="c", free=FALSE, values=0 )
AEFit     <- mxRun(twinAE)
AESumm   <- summary(AEFit)

# Generate AE Model Output
estVA_AE     <- mxEval(a*a, AEFit)               # additive genetic variance, a^2
estVE_AE     <- mxEval(e*e, AEFit)               # unique environmental variance, e^2
estVP_AE     <- (estVA_AE+estVE_AE)                    # total variance
estPropVA_AE <- estVA_AE/estVP_AE                      # standardized additive genetic variance
estPropVE_AE <- estVE_AE/estVP_AE                     # standardized unique environmental variance
estAE_AE     <- rbind(cbind(estVA_AE,estVE_AE),        # table of estimates
                   cbind(estPropVA_AE,estPropVE_AE))
LL_AE     <- mxEval(objective, AEFit)         # likelihood of AE model

#DF_AE = msize+dsize-nrow(AEFit@output$standardErrors)

#mychi_AE = LL_AE - LL_ACE
#mychi_DF_AE = DF_AE - DF_ACE
#mychi_p_AE = 1 - pchisq(mychi_AE, mychi_DF_AE)

#p_AE[edge] <- mychi_p_AE

#Generate CE Model
twinCE   <- twinACE
twinCE   <- omxSetParameters(twinCE, labels="a", free=FALSE, values=0 )
CEFit     <- mxRun(twinCE)
CESumm   <- summary(CEFit)

#Generate CE Model Output
estVC_CE     <- mxEval(c*c, CEFit)               # shared environmental variance, c^2
estVE_CE     <- mxEval(e*e, CEFit)               # unique environmental variance, e^2
estVP_CE     <- (estVC_CE+estVE_CE)                    # total variance
estPropVC_CE <- estVC_CE/estVP_CE                      # standardized additive genetic variance
estPropVE_CE <- estVE_CE/estVP_CE                     # standardized unique environmental variance
estAE_CE     <- rbind(cbind(estVC_CE,estVE_CE),        # table of estimates
                      cbind(estPropVC_CE,estPropVE_CE))
LL_CE     <- mxEval(objective, CEFit)         # likelihood of AE model

DF_CE = msize+dsize-nrow(CEFit@output$standardErrors)
mychi_CE = LL_CE - LL_ACE
mychi_DF_CE = DF_CE - DF_ACE
mychi_p_CE = 1 - pchisq(mychi_CE, mychi_DF_CE)
p_CE[edge] <- mychi_p_CE

#Generate AC Model
#twinAC   <- twinACE
#twinAC   <- omxSetParameters(twinAC, labels="e", free=FALSE, values=0 )
#ACFit     <- mxRun(twinAC)
#ACSumm   <- summary(ACFit)

#Generate AC Model Output
#estVA_AC     <- mxEval(a*a, ACFit)               # additive genetic variance, a^2
#estVC_AC     <- mxEval(c*c, ACFit)               # shared environmental variance, c^2
#estVP_AC     <- (estVA_AC+estVC_AC)                    # total variance
#estPropVA_AC <- estVA_AC/estVP_AC                      # standardized additive genetic variance
#estPropVC_AC <- estVC_AC/estVP_AC                     # standardized unique environmental variance
#estAC_AC     <- rbind(cbind(estVA_AC,estVC_AC),        # table of estimates
#                      cbind(estPropVA_AC,estPropVC_AC))
#LL_AC     <- mxEval(objective, ACFit)         # likelihood of AE model

#DF_AC = msize+dsize-nrow(ACFit@output$standardErrors)

#mychi_AC = LL_AC - LL_ACE
#mychi_DF_AC = DF_AC - DF_ACE
#mychi_p_AC = 1 - pchisq(mychi_AC, mychi_DF_AC)

#p_AC[edge] <- mychi_p_AC

options('digits' = 5)
compValues = mxCompare(twinACEFit, c(CEFit, AEFit))
Pvals[edge,] = compValues$p[2:3]

}
fileNameSave = sprintf("%s_%s_%s.txt", parcellation, tract, weights)
write.table(data.frame(heritabilityA[1:10]), fileNameSave, sep="\t")
