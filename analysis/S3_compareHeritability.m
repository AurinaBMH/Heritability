function S3_compareHeritability(parcellation,tractography,weight2,densThreshold,cvMeasure)

if strcmp(parcellation, 'HCP') || strcmp(parcellation, 'cust100')
    numSubc = 10;
elseif strcmp(parcellation, 'aparcaseg')
    numSubc = 7;
elseif strcmp(parcellation, 'cust250')
    numSubc = 15;
end
    
load(sprintf('twinEdges_%s_%s_%s_%s%d.mat', parcellation, tractography, weight2, cvMeasure, round(densThreshold*100))); 
heritFile = sprintf('heritabilityACE_%s_%s_%s_%s%d.txt', parcellation, tractography, weight2, cvMeasure, round(densThreshold*100)); 
heritabilityACE = importHeritabilityResult(heritFile); 

numNodes = size(groupAdjlog,1); 
% reshape heritability vector into the matrix for connected edges get indexes on diagonal
heritMatrix = zeros(numNodes,numNodes);
% mask upper half to get indexes of existing links
C = maskuHalf(groupAdjlog); 
% combine values to a vector for reshaping
% in heritability variable 1st column A, 2nd column C, 3rd column E
heritMatrix(C==1) = heritabilityACE.heritabilityA; % assign heritability values to those edges
% make a full matrix
heritMatrix = heritMatrix+heritMatrix'; 
nodeData = degrees_und(groupAdjlog); 
% make a curve plot for the whole brain
RichClubHuman(groupAdjlog,heritMatrix, nodeData); 
title('Whole brain')
ylabel('Mean edge heritability')
ylim([0.4 0.65])

% get values for the left cortex as in CGE analysis
% get values for the left and rifht cortex
LC = 1:numNodes/2-numSubc; 

groupAdjlogLC = groupAdjlog(LC,LC);
heritMatrixLC = heritMatrix(LC,LC);
nodeDataLC = nodeData(LC); 
RichClubHuman(groupAdjlogLC,heritMatrixLC, nodeDataLC)
title('Left cortex')
ylabel('Mean edge heritability')
ylim([0.4 0.65])

% get values for the left and rifht cortex
LC = 1:numNodes/2-numSubc; 
RC = numNodes/2+1:numNodes-numSubc; 
LCRC = [LC,RC]; 

groupAdjlogLRC = groupAdjlog(LCRC,LCRC);
heritMatrixLRC = heritMatrix(LCRC,LCRC);
nodeDataLRC = nodeData(LCRC);
RichClubHuman(groupAdjlogLRC,heritMatrixLRC, nodeDataLRC)
title('Cortex')
ylabel('Mean edge heritability')
ylim([0.4 0.65])
end

