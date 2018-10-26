% This script is to be used for heritability analyses 
%---------------------------------------------------------------------
% 1. First step is to create covariates that include MZ twins, DZ twins and
% their siblings. S0_giveCovariatesDWI does that; This steps should be done
% only ince, uncomment the following line if you want to run it. 
% S0_giveCovariatesDWI 
%---------------------------------------------------------------------
%---------------------------------------------------------------------
% 2. Second step is to create the variable that contains edges for the
% heritability analysis: For this a group connectome is generated using 
% streamline count as weight and using this mask edges for each subject are
% selected using FA as weight given it's more biological nature. 
% This is one using function S1_giveEdges with defined options.
%---------------------------------------------------------------------
% run S1_giveEdges
% choose options for the data to be loaded
whatDATA = 'HCP'; % GenCog or HCP data
parcellation = 'HCP'; % type of parcellation
tractography = 'iFOD2'; % tractography algorithm - FACT vs iFOD2
weight1 = 'standard'; % weight type fir the initial group connectome generation
brainPart = 'wholeBrain'; % what part of the brain to consider 
strRem1 = 10; % connections with the <strRem1 streamlines will be excluded from the connectome generation; 
densThreshold = 0.15; % density of the group connectome
groupConn = 'CVmeasure'; % type of the group connectome
cvMeasure = 'strength'; % type of measure to quantify the the consistency of the weights.  
consThr = 0.6; % connections with the consistency <consThr will not be considered for the group connectome

weight2 = 'FA'; % type of weight that will be used in the heritability analysis
strRem2 = 0; % no connecions to be set to 0. 

S1_giveEdges(whatDATA,parcellation,tractography,brainPart,weight1,strRem1,densThreshold, groupConn,cvMeasure,consThr,weight2,strRem2)

%---------------------------------------------------------------------
% 3. Third step is to run heritabiloty analysis in R using the same options
% - at this point you need to open Rstudio and run the script from there
%---------------------------------------------------------------------

S2_ACEheritability(parcellation, tractography, weight2, cvMeasure, densThreshold)
%[RESULT, STATUS, MSG] = evalR('source('~/R functions/makerm function.R')');  
system('"/Library/Frameworks/R.framework/Resources/Rscript" /Users/Aurina/GoogleDrive/Genetics_connectome/Heritability/code/analysis/S2_ACEheritability(parcellation, tractography, weight2, cvMeasure, densThreshold)')
%---------------------------------------------------------------------
% 4. Fourth step is to plot heritability as a function of 
%---------------------------------------------------------------------

S3_compareHeritability(parcellation,tractography,weight2,densThreshold,cvMeasure)





