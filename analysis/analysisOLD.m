
close all;
clear all;

% choose options for the data to be loaded
whatDWI = 'HCP';
parc = 'HCP';
tract = 'iFOD2';
weight1 = 'standard';
brainPart = 'wholeBrain'; 
strRem = 10;

[A, matrices, coordinates, avWeight, SUBjects] = giveConnDATA(whatDWI,parc,tract,weight1,brainPart,strRem); 
load('twinCovariatesDWI.mat')

subjects = vertcat(MZ_ID(:), DZ_ID(:)); 
subjects(isnan(subjects)) = []; 

[SUBS, subIND] = intersect(SUBjects, subjects);

coordinatesMask = coordinates(subIND);
connectomesMask = matrices(subIND);
%SUBS = ConnMask.SUBS(subIND);


% make a group matrix and select only edges existing in the group matrix
numNodes = size(coordinatesMask{1},1);
numSubj = size(coordinatesMask,2);

%-----------------------------------------------------------------
% Make a group matrix using streamline count as weight
%-----------------------------------------------------------------
groupConn = 'CVmeasure'; %'CVmeasure'; %  lengthCV, consistency
cvMeasure = 'strength';
giveRC=false; 
densThreshold = 0.15; 
consThr = 0.6; 

[Gr] = giveMeRichClub(matrices, coordinates, groupConn ,densThreshold, giveRC, cvMeasure, consThr);
groupAdjlog = logical(Gr); 
dens = density_und(Gr); 

%-----------------------------------------------------------------
% Select edges using FA as weight
%-----------------------------------------------------------------
weight2 = 'FA';
strRem = 0; 
[A, matrices, coordinates, avWeight, SUBjects] = giveConnDATA(whatDWI,parc,tract,weight2,brainPart,strRem); 
coordinatesMask = coordinates(subIND);
connectomesMask = matrices(subIND);

% get indeces for non-zero elements
m = maskuHalf(connectomesMask{1}.*groupAdjlog);
m(groupAdjlog==0) = NaN;
Wind = find(~isnan(m)); 

edgeMatr = zeros(sum(groupAdjlog(:))/2, numSubj);
for s=1:numSubj
    m = maskuHalf(connectomesMask{s}.*groupAdjlog);
    m(groupAdjlog==0) = NaN;
    edgeMatr(:,s) = m(~isnan(m));
end

numEdges = size(edgeMatr,1);
Output_MZ = zeros(size(MZ_ID,1),size(MZ_ID,2),numEdges);

for i=1:size(MZ_ID,1)
    ind1 = find(SUBS==MZ_ID(i,1));
    ind2 = find(SUBS==MZ_ID(i,2));
    ind3 = find(SUBS==MZ_ID(i,3));
    ind4 = find(SUBS==MZ_ID(i,4));
    for edg = 1:numEdges
        Output_MZ(i,1,edg) = edgeMatr(edg,ind1);
        Output_MZ(i,2,edg) = edgeMatr(edg,ind2);
        
        if ~isempty(ind3)
            Output_MZ(i,3,edg) = edgeMatr(edg,ind3);
        else
            Output_MZ(i,3,edg) = NaN; 
        end
        
        if ~isempty(ind4)
            Output_MZ(i,4,edg) = edgeMatr(edg,ind4);
        else
            Output_MZ(i,4,edg) = NaN; 
        end
    end
end

Output_DZ = zeros(size(DZ_ID,1),size(DZ_ID,2), numEdges);
for i=1:size(DZ_ID,1)
    ind1 = find(SUBS==DZ_ID(i,1));
    ind2 = find(SUBS==DZ_ID(i,2));
    ind3 = find(SUBS==DZ_ID(i,3));
    ind4 = find(SUBS==DZ_ID(i,4));
    for edg = 1:numEdges
        Output_DZ(i,1,edg) = edgeMatr(edg,ind1);
        Output_DZ(i,2,edg) = edgeMatr(edg,ind2);
        
        if ~isempty(ind3)
            Output_DZ(i,3,edg) = edgeMatr(edg,ind3);
        else
            Output_DZ(i,3,edg) = NaN; 
        end
        
        if ~isempty(ind4)
            Output_DZ(i,4,edg) = edgeMatr(edg,ind4);
        else
            Output_DZ(i,4,edg) = NaN; 
        end
    end
end

cd ..
cd ('output')
fileName = sprintf('twinEdges_%s_%s_%sTEST%d.mat', parcellation, tract, weight, dens*100);
save(fileName, 'Output_MZ', 'Output_DZ');

% count how many pairs of subject we'd need to exclude if a pair is
% excluded if any twin has a value of 0 on the given edge
edgeDZ = zeros(numEdges,1);
edgeMZ = zeros(numEdges,1);
for ed=1:numEdges
    t1DZ = find(Output_DZ(:,1,ed)==0);
    t2DZ = find(Output_DZ(:,2,ed)==0);
    
    t1MZ = find(Output_MZ(:,1,ed)==0);
    t2MZ = find(Output_MZ(:,2,ed)==0);
    
    if ~isempty(t1DZ) || ~isempty(t2DZ)
        
        edgeDZ(ed) = max(length(t1DZ),length(t2DZ));
        
    end
    
    if ~isempty(t1MZ) || ~isempty(t2MZ)
        
        edgeMZ(ed) = max(length(t1MZ),length(t2MZ));
        
    end
end

figure;
subplot(1,2,1); histogram(edgeDZ); title('DZ'); ylabel('Number of edges'); xlabel('Number of twin pairs to exclude')
subplot(1,2,2); histogram(edgeMZ);title('MZ'); ylabel('Number of edges'); xlabel('Number of twin pairs to exclude')

% label edges based on hub-ness
deg = degrees_und(groupAdjlog);
k=1; 
for khub=10:10:180

isHub = deg>khub;

% label all edges between nodes as rich, feeder or peripheral
mask = zeros(numNodes, numNodes);
mask(isHub, isHub) = 3;
mask(isHub, ~isHub) = 2;
mask(~isHub, isHub) = 2;
mask(~isHub, ~isHub) = 1;

%selectMask = mask.*groupAdjlog;

selectMask = maskuHalf(mask.*groupAdjlog);
selectMask(groupAdjlog==0) = NaN;
edgeLabel = selectMask(~isnan(selectMask));

% load h^2 data from heritability analyses and plot the distributions
% for those three groups;

heritability = importfile('heritabilityACE.txt'); 
dataCell = cell(3,1);
dataCell{1} = heritability(edgeLabel==3,1); % rich;
dataCell{2} = heritability(edgeLabel==2,1); % feeder
dataCell{3} = heritability(edgeLabel==1,1); % peripheral;

JitteredParallelScatter(dataCell); xticks([1 2 3]); 
xticklabels({sprintf('rich %d', length(dataCell{1})), ...
   sprintf('feeder %d',length(dataCell{2})), sprintf('peripheral %d', length(dataCell{3}))})
title(sprintf('hub threshold %d', khub))


[pRF(k),~,statsRF.tstat.zval.zval] = ranksum(dataCell{1}, dataCell{2}); % rich VS feeder
tRF(k) = statsRF.tstat.zval.zval; 
[pRP(k),~,statsRP.tstat.zval.zval] = ranksum(dataCell{1}, dataCell{3}); % rich VS peripheral
tRP(k) = statsRP.tstat.zval.zval; 
[pFP(k),~,statsFP.tstat.zval.zval] = ranksum(dataCell{2}, dataCell{3}); % feeder vs peripheral
tFP(k) = statsFP.tstat.zval.zval; 
k=k+1; 
end
% reshape heritability vector into the matrix for connected edges
% get indexes on diagonal
K = zeros(numNodes,numNodes);
K(logical(eye(380))) = 1; 
Eind = find(K==1); 
% mask upper half to get indexes of existing and non-existing links
C = maskuHalf(groupAdjlog); 
Qind = find(C==0); 
%Wind = find(C==1); 

% combine values to a vector for reshaping
% in heritability variable 1st column A, 2nd column C, 3rd column E
W(Qind) = 0; W(Wind) = heritability(:,1); W(Eind) = 0; 
% reshape
heritMatrix = reshape(W,380,380); 
% make a full matrix
heritMatrix = heritMatrix+heritMatrix'; 

% make a curve plot
RichClubHuman(groupAdjlog,heritMatrix)



