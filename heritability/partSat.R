mylabels=c("twin1","twin2")
mxFitFunctionMultigroup(c("MZsat","DZsat"))
fitFunction <- mxFitFunctionML(rowDiagnostics=TRUE)

MZsat <- mxModel("MZsat",
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(0,0), name="expMeanMZ" ),
                 mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, .5, name="CholMZ" ),
                 mxAlgebra( CholMZ %*% t(CholMZ), name="expCovMZ"),
                 mxData( mzData_measure, type="raw"), mxFitFunctionML(),
                 mxExpectationNormal( "expCovMZ", "expMeanMZ", mylabels))

DZsat <- mxModel("DZsat",
                 mxMatrix( type = "Full", nrow=1, ncol=2, free=T, c(0,0), name="expMeanDZ" ),
                 mxMatrix( type = "Lower", nrow=2, ncol=2, free=T, .5, name="CholDZ" ),
                 mxAlgebra( CholDZ %*% t(CholDZ), name="expCovDZ"),
                 mxData( dzData_measure, type="raw"), mxFitFunctionML(),
                 mxExpectationNormal( "expCovDZ", "expMeanDZ", mylabels))

mytwinSatModel <- mxModel(model='twinSat', MZsat, DZsat, mxFitFunctionMultigroup(c('MZsat', 'DZsat')))
saturatedOut <- mxRun(mytwinSatModel)
summary(saturatedOut)
