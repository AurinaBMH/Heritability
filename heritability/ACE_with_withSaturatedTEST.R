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