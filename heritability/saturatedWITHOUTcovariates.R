# example

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

data_covar$DZ.age = data_covar$DZ.age/100
data_covar$MZ.age = data_covar$MZ.age/100

mxOption(NULL,"Default optimizer","SLSQP")
myDataMZ<-data.frame(data$Output.MZ)
myDataDZ<-data.frame(data$Output.DZ)

colnames(myDataMZ) <- c('twin1', 'twin2')
colnames(myDataDZ) <- c('twin1', 'twin2')
selVars <- c('twin1','twin2')

summary(myDataMZ)
summary(myDataDZ)
colMeans(myDataMZ,na.rm=TRUE)  
#na.rm means ignore NA values (non-numeric)
colMeans(myDataDZ,na.rm=TRUE)
cov(myDataMZ,use="complete")
# "complete" specifies use only cases with data in all columns
cov(myDataDZ,use="complete")

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
                                  mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(0,0), name="expMeanMZ" ),
                                  mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, .5, name="CholMZ" ),
                                  mxAlgebra( CholMZ %*% t(CholMZ), name="expCovMZ"),
                                  mxData( myDataMZ, type="raw"), mxFitFunctionML(),
                                  mxExpectationNormal( "expCovMZ", "expMeanMZ", mylabels))
                 
DZsat <- mxModel("DZsat",
                                  mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(0,0), name="expMeanDZ" ),
                                  mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, .5, name="CholDZ" ),
                                  mxAlgebra( CholDZ %*% t(CholDZ), name="expCovDZ"),
                                  mxData( myDataDZ, type="raw"), mxFitFunctionML(),
                                  mxExpectationNormal( "expCovDZ", "expMeanDZ", mylabels))

#---------------------------------------------------------------------------------------------
#            Saturated_twin_model
#            DVM Bishop, 12th March 2010, based on OpenMxUsersGuide, p. 16
#---------------------------------------------------------------------------------------------

colnames(myDataMZ)=c("twin1","twin2")
colnames(myDataDZ)=colnames(myDataMZ)
mylabels=c("twin1","twin2")

# Model specification starts here
SatModel <- mxModel("twinSat", MZsat, DZsat, mxFitFunctionMultigroup(c('MZsat', 'DZsat')))
# adds together likelihoods for MZ and DZ groups
# evaluate expression from mxAlgebra, i.e. both submodels together
#---------------------------------------------------------------------------------------------------------------------------
SatModelFit <- mxTryHard(SatModel) #The mxRun command evaluates the model.
summary(SatModelFit)
#-------------------------------------------------------------------------------------------------------------------
 
# at the moment doesn't run when covariates are included in the saturated model
# (2) Specify and Run Sub-Model1: equating means across twin order within zygosity group
# note: Within zyg group we equate Tw1 and Tw2 Means by using one label iso two
# ------------------------------------------------------------------------------------------------------------
Sub1Model		<- SatModel
Sub1Model$MZsat$expMeanMZ	<- mxMatrix(type="Full", nrow=1, ncol=2, free=T, c(0,0), label="mMZ", name="expMeanMZ")
Sub1Model$DZsat$expMeanDZ	<- mxMatrix(type="Full", nrow=1, ncol=2, free=T, c(0,0), label="mDZ", name="expMeanDZ")

Sub1Fit			<- mxTryHard(Sub1Model)
(Sub1Sum		<- summary(Sub1Fit))

pVals = numeric(4)
# Test birth order effect on the means
comparison01 <- mxCompare(SatModelFit, Sub1Fit) 
pVals[1] <- comparison01$p[2]


# (3) Specify and Run Sub-Model2: equating means across twin order AND zygosity group
# note: From previous model we equate all means by using one label across zyg groups 
# ----------------------------------------------------------------------------------------------------------

Sub2Model			<- Sub1Model
Sub2Model$MZsat$expMeanMZ	<- mxMatrix(type="Full", nrow=1, ncol=2, free=T, c(0,0), label="mean", name="expMeanMZ")
Sub2Model$DZsat$expMeanDZ	<- mxMatrix(type="Full", nrow=1, ncol=2, free=T, c(0,0), label="mean", name="expMeanDZ")
Sub2Fit 			<- mxTryHard(Sub2Model)
(Sub2Sum		 	<- summary(Sub2Fit))

# Test zygosity effect on the means
comparison12 <-mxCompare(Sub1Fit, Sub2Fit)  
pVals[2] <- comparison12$p[2]

#***********************************************************
# Now we will do the same equality testing for the variances
#***********************************************************

# (4) Specify and Run Sub-Model3: equating Variances across twin order in MZ and DZ group
# note: From previous model we equate Var of Tw 1 and tw by by using one iso two labels
# ----------------------------------------------------------------------------------------------------------
nv<- 1

Sub3Model			<- Sub2Model
Sub3Model$MZsat$VarMZtw1	<- mxMatrix(type="Full", nrow=nv, ncol=nv, free=TRUE, c(0,0), label=c("VMZ"),name="VarMZtw1")
Sub3Model$MZsat$VarMZtw2	<- mxMatrix(type="Full", nrow=nv, ncol=nv, free=TRUE, c(0,0), label=c("VMZ"),name="VarMZtw2")
Sub3Model$DZsat$VarDZtw1	<- mxMatrix(type="Full", nrow=nv, ncol=nv, free=TRUE, c(0,0), label=c("VDZ"),name="VarDZtw1")
Sub3Model$DZsat$VarDZtw2	<- mxMatrix(type="Full", nrow=nv, ncol=nv, free=TRUE, c(0,0), label=c("VDZ"),name="VarDZtw2")
Sub3Fit			<- mxTryHard(Sub3Model)
(Sub3Sum		<- summary(Sub3Fit))

#Test birth order effect on the Variance
comparison32<-mxCompare(Sub3Fit,Sub2Fit) 
pVals[3] <- comparison32$p[2]

# (5) Specify and Run Sub-Model4: equating Variances across twin order AND zygosity group
# note: From previous model we equate all variances by using one label across zyg groups 
# ---------------------------------------------------------------------------------------------------------------

Sub4Model			<- Sub3Model
Sub4Model$MZsat$VarMZtw1	<- mxMatrix(type="Full", nrow=nv, ncol=nv, free=TRUE, c(0,0), label=c("V"),name="VarMZtw1")
Sub4Model$MZsat$VarMZtw2	<- mxMatrix(type="Full", nrow=nv, ncol=nv, free=TRUE, c(0,0), label=c("V"),name="VarMZtw2")
Sub4Model$DZsat$VarDZtw1	<- mxMatrix(type="Full", nrow=nv, ncol=nv, free=TRUE, c(0,0), label=c("V"),name="VarDZtw1")
Sub4Model$DZsat$VarDZtw2	<- mxMatrix(type="Full", nrow=nv, ncol=nv, free=TRUE, c(0,0), label=c("V"),name="VarDZtw2")
Sub4Fit			<- mxTryHard(Sub4Model)
(Sub4Sum		<- summary(Sub4Fit))

#Test difference in Variance across zygosity group
comparison34<- mxCompare(Sub3Fit,Sub4Fit) 
pVals[4] <- comparison34$p[2]
