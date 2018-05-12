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

nv	<- 1		# number of variables for a twin = 1 in Univariate
ntv	<- 2*nv	# number of variables for a pair = 2* 1 for Univariate


options(warn=1)
load('dataTEST.RData')
# all edges

mzDataALL<-data.frame(data$Output.MZ, data_covar$MZ.age[,1], data_covar$MZ.age[,2], data_covar$MZ.sex[,1], data_covar$MZ.sex[,2])
dzDataALL<-data.frame(data$Output.DZ, data_covar$DZ.age[,1], data_covar$DZ.age[,2], data_covar$DZ.sex[,1], data_covar$DZ.sex[,2])


colnames(mzData) <- c('twin1', 'twin2', 'ageT1MZ', 'ageT2MZ', 'sexT1MZ', 'sexT2MZ')
colnames(dzData) <- c('twin1', 'twin2', 'ageT1DZ', 'ageT2DZ', 'sexT1DZ', 'sexT2DZ')

names (mzDataALL)
str(mzDataALL)
summary(mzDataALL)
describe(mzDataALL)

colnames(dzData) <- c('iq1','iq2')
colnames(mzData) <- c('iq1','iq2')

mzData		<- data$Output.MZ
dzData		<- data$Output.DZ

describe(mzData)
describe(dzData)

colMeans(mzData,na.rm=TRUE)
cov(mzData,use="complete")
colMeans(dzData,na.rm=TRUE)
cov(dzData,use="complete")


# -----------------------------------------------------------------------
# (1) Specify and Run Saturated Model (Cholesky Decomposition)
# -----------------------------------------------------------------------
# Specify Matrices
CholMZ	<-mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=T, values=5, name="lowMZ" )
CholDZ	<-mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=T, values=5, name="lowDZ" )
MZCov	<-mxAlgebra( expression=lowMZ %*% t(lowMZ), name="expCovMZ" )
DZCov	<-mxAlgebra( expression=lowDZ %*% t(lowDZ), name="expCovDZ" )

MZMeans	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=90, labels=c("Mmz1","Mmz2"), name="expMeanMZ" )
DZMeans	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=90, labels=c("Mdz1","Mdz2"), name="expMeanDZ" )

# Algebra's needed for picking out elements from Cov models (used later for equality constraints)
MZexpVar	<-mxAlgebra( expression= t(diag2vec(expCovMZ)), name="expVarMZ")
MZexpVartw1	<-mxAlgebra( expression= expVarMZ[1,1], 	name="expVarMZt1")
MZexpVartw2	<-mxAlgebra( expression= expVarMZ[1,2], 	name="expVarMZt2")
DZexpVar	<-mxAlgebra( expression= t(diag2vec(expCovDZ)), name="expVarDZ")
DZexpVartw1	<-mxAlgebra( expression= expVarDZ[1,1], 	name="expVarDZt1")
DZexpVartw2	<-mxAlgebra( expression= expVarDZ[1,2], 	name="expVarDZt2")

# Algebra's needed for standardizing the covariances
matI	<-mxMatrix( type="Iden", nrow=ntv, ncol=ntv, name="I")
MZcor	<-mxAlgebra( expression= solve(sqrt(I*expCovMZ)) %*% expCovMZ %*% solve(sqrt(I*expCovMZ)), name="rMZ" )
DZcor	<-mxAlgebra( expression= solve(sqrt(I*expCovDZ)) %*% expCovDZ %*% solve(sqrt(I*expCovDZ)), name="rDZ" )

# Data objects for Multiple Groups
dataMZ	<-mxData(mzData, type="raw")
dataDZ	<-mxData(dzData, type="raw")

# Objective objects for Multiple Groups
# mxExpectationNormal: Objective functions which uses Full?Information maximum
# likelihood, the preferred method for raw data.
# Objective functions are functions for which free parameter values are chosen
# such that the value of the objective function is minimized
objMZ	<-mxExpectationNormal( covariance="expCovMZ", means="expMeanMZ", dimnames=selVars )
objDZ	<-mxExpectationNormal( covariance="expCovDZ", means="expMeanDZ", dimnames=selVars )

fitFunction <- mxFitFunctionML()

# Combine Groups
groupMZ	<-mxModel("MZ", CholMZ, MZCov, MZMeans, MZexpVar, MZexpVartw1,MZexpVartw2, matI, MZcor, dataMZ, objMZ, fitFunction )
groupDZ	<-mxModel("DZ", CholDZ, DZCov, DZMeans, DZexpVar, DZexpVartw1,DZexpVartw2, matI, DZcor, dataDZ, objDZ, fitFunction )
minus2ll  	<-mxAlgebra( MZ.objective + DZ.objective, name="minus2LL" )
obj		<-mxFitFunctionAlgebra("minus2LL")
SatModel   	<-mxModel( "Sat", minus2ll, obj, groupMZ, groupDZ )

SatFit	<- mxRun(SatModel)
(SatSum 	<- summary(SatFit))

# Generate some output
SatFit$MZ$expCovMZ
SatFit$DZ$expCovDZ
SatFit$MZ$rMZ
SatFit$DZ$rDZ

mxEval(MZ.expMeanMZ, SatFit)
mxEval(MZ.expCovMZ, SatFit)
mxEval(DZ.expMeanDZ, SatFit)
mxEval(DZ.expCovDZ, SatFit)




















# load data - data is not attached to the question due to file format differences
# -----------------------------------------------------------------------
options(warn=1)
load('dataTEST.RData')
# all edges
numEdges = 1

MeasureMZ = 'Edge weight MZ'
MeasureDZ = 'Edge weight DZ'

mzData_measure<-data$Output.MZ
dzData_measure<-data$Output.DZ

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

pictureName = "testPlot.png"
png (pictureName)
split.screen(c(1,2))
screen(1)
plot(mzData_measure,main=MeasureMZ,xlim=c(min(mzData_measure, dzData_measure), max(mzData_measure, dzData_measure)), ylim=c(min(mzData_measure, dzData_measure), max(mzData_measure, dzData_measure)))
screen(2)
plot(dzData_measure,main=MeasureDZ,xlim=c(min(mzData_measure, dzData_measure), max(mzData_measure, dzData_measure)), ylim=c(min(mzData_measure, dzData_measure), max(mzData_measure, dzData_measure)))
dev.off() # to complete the writing process and return output to your monitor
# -----------------------------------------
mzData<-data.frame(data$Output.MZ, data_covar$MZ.age[,1], data_covar$MZ.age[,2], data_covar$MZ.sex[,1], data_covar$MZ.sex[,2])
dzData<-data.frame(data$Output.DZ, data_covar$DZ.age[,1], data_covar$DZ.age[,2], data_covar$DZ.sex[,1], data_covar$DZ.sex[,2])

colnames(mzData) <- c('twin1', 'twin2', 'ageT1MZ', 'ageT2MZ', 'sexT1MZ', 'sexT2MZ')
colnames(dzData) <- c('twin1', 'twin2', 'ageT1DZ', 'ageT2DZ', 'sexT1DZ', 'sexT2DZ')
selVars <- c('twin1','twin2')

# Fit saturated model where relationships between twin pairs are ignored - models for MZ and DZ twins are identical.
#--------------------------------------------------------------------------------------------
# Model specification starts here
mylabels=c("twin1","twin2")
mxFitFunctionMultigroup(c("MZsat","DZsat"))
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

saturatedOut <- mxRun(mytwinSatModel, intervals = T)
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
twinACEFit <- mxRun(twinACE)
# running the same model twice - increases heritability for the 2nd edge
twinACEFit <- mxRun(twinACEFit)

summary(twinACEFit)

# compare to saturated model
compValues = mxCompare(saturatedOut, twinACEFit)
