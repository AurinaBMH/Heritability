
#----------------------------------------------------------------------------------------------------------------------
#                                                             Twin simulate
# DVM Bishop, 11th March 2010, Based on script in OpenMXUserGuide, p 15
#----------------------------------------------------------------------------------------------------------------------
require(OpenMx)   # not needed for the simulation, but will be needed when we come to model specification
require(MASS)    # needed for multivariate random number generation
library(R.matlab)
set.seed(200)        # specified seed ensures same random number set generated on each run

setwd("~/GoogleDrive/Genetics_connectome/Heritability/data/output")
data_covar = readMat("~/GoogleDrive/Genetics_connectome/Heritability/data/general/twinCovariatesDWI.mat")
data = readMat("~/GoogleDrive/Genetics_connectome/Heritability/data/output/twinDataTEST.mat")

mxOption(NULL,"Default optimizer","SLSQP")
data_covar$DZ.age = data_covar$DZ.age/100
data_covar$MZ.age = data_covar$MZ.age/100

myDataMZ<-data.frame(data$Output.MZ, data_covar$MZ.age[,1], data_covar$MZ.age[,2], data_covar$MZ.sex[,1], data_covar$MZ.sex[,2])
myDataDZ<-data.frame(data$Output.DZ, data_covar$DZ.age[,1], data_covar$DZ.age[,2], data_covar$DZ.sex[,1], data_covar$DZ.sex[,2])

colnames(myDataMZ) <- c('twin1', 'twin2','ageT1MZ', 'ageT2MZ', 'sexT1MZ', 'sexT2MZ')
colnames(myDataDZ) <- c('twin1', 'twin2','ageT1DZ', 'ageT2DZ', 'sexT1DZ', 'sexT2DZ')
selVars <- c('twin1','twin2')



summary(myDataMZ)
summary(myDataDZ)
colMeans(myDataMZ,na.rm=TRUE)  
#na.rm means ignore NA values (non-numeric)
colMeans(myDataDZ,na.rm=TRUE)
cov(myDataMZ,use="complete")
# "complete" specifies use only cases with data in all columns
cov(myDataDZ,use="complete")

cMZ = cov(data$Output.MZ,use="complete")
cMZ[1,2] = 0
# "complete" specifies use only cases with data in all columns
cDZ = cov(data$Output.DZ,use="complete")
cDZ[1,2] = 0

# do scatterplots for MZ and DZ
split.screen(c(1,2))        # split display into two screens side by side
screen(1)
plot(data$Output.DZ,main='DZ')    # main specifies overall plot title
screen(2)
plot(data$Output.MZ,main='MZ')
#use drag and drop to resize the plot window if necessary

#--------------------------------------------------------------------------------------------
#  DZ submodel, illustrating compact formatting
#--------------------------------------------------------------------------------------------
mylabels=c("twin1","twin2")

MZsat <- mxModel("MZsat",
                 
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(0,0), labels =c("b0_mz1","b0_mz2"), name="Intercepts" ),
                 mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                 mxMatrix( type="Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ"), name="MZDefVars"),
                 mxAlgebra( expression=Intercepts + beta %*% MZDefVars, name="expMeanMZ"),
                 mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, 0.5, name="CholMZ" ),
                 mxAlgebra( CholMZ %*% t(CholMZ), name="expCovMZ"),
                 mxData( myDataMZ, type="raw"), mxFitFunctionML(),
                 mxExpectationNormal( "expCovMZ", "expMeanMZ", mylabels))

DZsat <- mxModel("DZsat",
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(0,0), labels=c("b0_dz1","b0_dz2"), name="Intercepts" ),
                 mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                 mxMatrix( type="Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ"), name="DZDefVars"),
                 mxAlgebra( expression=Intercepts + beta %*% DZDefVars, name="expMeanDZ"),
                 mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, 0.5, name="CholDZ" ),
                 mxAlgebra( CholDZ %*% t(CholDZ),name="expCovDZ"),
                 mxData( myDataDZ, type="raw"), mxFitFunctionML(),
                 mxExpectationNormal( "expCovDZ", "expMeanDZ", mylabels))

mxMatrix( type="Full", nrow=2, ncol=2, free=F, label=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ"), name="MZDefVars"),
mxAlgebra( expression=twinACE.expMean + twinACE.beta %*% MZDefVars, name="expMeanMZ"),
mxExpectationNormal(covariance="twinACE.expCovMZ", means="expMeanMZ", dimnames=selVars ),
mxFitFunctionML()),


#---------------------------------------------------------------------------------------------
#            Saturated_twin_model
#            DVM Bishop, 12th March 2010, based on OpenMxUsersGuide, p. 16
#---------------------------------------------------------------------------------------------

#colnames(myDataMZ)=c("twin1","twin2")
#colnames(myDataDZ)=colnames(myDataMZ)
mylabels=c("twin1","twin2")

# Model specification starts here
SatModel <- mxModel("twinSat", MZsat, DZsat, mxFitFunctionMultigroup(c('MZsat', 'DZsat')))
# adds together likelihoods for MZ and DZ groups
# evaluate expression from mxAlgebra, i.e. both submodels together
#---------------------------------------------------------------------------------------------------------------------------
SatModelFit <- mxRun(SatModel) #The mxRun command evaluates the model.
summary(SatModelFit)

#Equate intercepts across twin order:
sub1 <- omxSetParameters(model=SatModel,labels=c("b0_mz1","b0_mz2"),newlabels="b0_mz",name="Submodel1")
sub1 <- omxSetParameters(model=sub1,labels=c("b0_dz1","b0_dz2"),newlabels="b0_dz",name="Submodel1")
sub1Fit <- mxRun(sub1)
mxCompare(SatModelFit, sub1Fit) # compare models

#Equate intercepts across zygosity:
sub2 <- omxSetParameters(model=sub1,labels=c("b0_mz","b0_dz"),newlabels="b0",name="Submodel2")
sub2Fit <- mxRun(sub2)
mxCompare(sub1Fit, sub2Fit) # compare models

#Equate variance across twin order:
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























                
sub3 <- omxSetParameters(model=sub2,labels=c("mzv1","mzv2"),newlabels="mzv",name="Submodel3")
sub3 <- omxSetParameters(model=sub3,labels=c("dzv1","dzv2"),newlabels="dzv",name="Submodel3")
sub3Fit <- mxRun(sub3)
mxCompare(sub2Fit, sub3Fit) # compare models

#Equate variance across zygosity:
sub4 <- omxSetParameters(model=sub3,labels=c("mzv","dzv"),newlabels="v",name="Submodel4")
sub4Fit <- mxRun(sub4)
mxCompare(sub3Fit, sub4Fit) # compare models
