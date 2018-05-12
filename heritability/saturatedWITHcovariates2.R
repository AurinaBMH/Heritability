
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
cMZ = cov(data$Output.MZ,use="complete")
# "complete" specifies use only cases with data in all columns
cDZ = cov(data$Output.DZ,use="complete")
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
                 mxMatrix(type="Symm",nrow=2,free=T,values=cMZ,labels=c("mzv1","mz","mzv2"),name="expCovMZ"),
                 mxData( myDataMZ, type="raw"), mxFitFunctionML(),
                 mxExpectationNormal( "expCovMZ", "expMeanMZ", mylabels))

DZsat <- mxModel("DZsat",
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(0,0), labels=c("b0_dz1","b0_dz2"), name="Intercepts" ),
                 mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                 mxMatrix( type="Full", nrow=2, ncol=2, free=F, labels=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ"), name="DZDefVars"),
                 mxAlgebra( expression=Intercepts + beta %*% DZDefVars, name="expMeanDZ"),
                 mxMatrix(type="Symm",nrow=2,free=T,values=cDZ,labels=c("dzv1","dz","dzv2"),name="expCovDZ"),
                 mxData( myDataDZ, type="raw"), mxFitFunctionML(),
                 mxExpectationNormal( "expCovDZ", "expMeanDZ", mylabels))

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