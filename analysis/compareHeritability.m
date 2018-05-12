% test the result

clear all; close all;
parcellation = 'HCPMMP1'; % 'HCPMMP1' , 'custom200';
tract = 'iFOD2';
tract2 = 'iFOD2';
sift = 'SIFT2';
groupConn = 'variance';
weight = 'FA'; % 'FA', 'standard'
dens = 0.05; 
onlyTwinAdj = true; 

heritability = importfile(sprintf('heritabilityACE_TESTHCP%d973.txt', dens*100)); 

cd ('data/general')

load('twinCovariatesDWI.mat')
cd ..
cd ('connectomes')


ConnMask = load(sprintf('%sANDfslatlas20_acpc_%s_%s_%s_structnets.mat', parcellation, tract, sift, weight));

twins = vertcat(MZ_ID(:), DZ_ID(:));
[~, twinIND] = intersect(ConnMask.SUBS, twins);
if onlyTwinAdj

coordinatesMask = ConnMask.COG(twinIND);
connectomesMask = ConnMask.ADJS(twinIND);
else
coordinatesMask = ConnMask.COG; %(twinIND);
connectomesMask = ConnMask.ADJS; %(twinIND);
end
    

% make a group matrix and select only edges existing in the group matrix
numNodes = size(coordinatesMask{1},1);
numSubj = size(coordinatesMask,2);

% make vectors for hemispheres for different parcellations based on the
% number of nodes (first half is always left, second half is always right)
hemiid = zeros(numNodes,1);
hemiid(1:numNodes/2) = 1;
hemiid(numNodes/2+1:numNodes) = 2;

dist = zeros(numNodes, numNodes, numSubj);
adjMatr = zeros(numNodes, numNodes, numSubj);

for s=1:numSubj
    dist(:,:,s) = pdist2(coordinatesMask{s}, coordinatesMask{s});
    %     connectomes = connectomesMask{s};
    %     connectomes(connectomes==0) = NaN;
    %     connectomesMask{s} = connectomes;
    adjMatr(:,:,s) = connectomesMask{s};
end
% take average of distance
avDist = mean(dist,3);

%-----------------------------------------------------------------
% Make three versions of group matrices - length, variance, consistency
%-----------------------------------------------------------------
if strcmp(groupConn, 'length')
    groupAdj_maskLength = fcn_group_average(adjMatr,avDist,hemiid);
    % replace zeros with NaNs;
    adjMatr(adjMatr==0) = NaN;
    % take average across subjectsclose all
    meanAdj = nanmean(adjMatr,3);
    groupAdj = meanAdj.*groupAdj_maskLength;
    groupAdj(isnan(groupAdj)) = 0;
elseif strcmp(groupConn, 'variance')
    [groupAdj, consist_var] = giveMeGroupAdj_variance(connectomesMask, dens);
elseif strcmp(groupConn, 'consistency')
    
    [groupAdj, groupDist, consist_cons] = giveMeGroupAdj_consistency(connectomesMask, distances,threshold);
end

groupAdjlog = logical(groupAdj);

%get indeces for existing values in group matrix (just on half so they can
%be used to ut heritability values back to the matrix)

m = maskuHalf(groupAdjlog);
m(groupAdjlog==0) = NaN;
Wind = find(~isnan(m)); 

% plot the distributions of heritability values based on hub-ness
% label edges based on hub-ness
% import heritability

deg = degrees_und(groupAdjlog);
k=1; 
for khub=20:10:100

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
K(logical(eye(numNodes))) = 1; 
Eind = find(K==1); 
% mask upper half to get indexes of existing and non-existing links
C = maskuHalf(groupAdjlog); 
Qind = find(C==0); 
%Wind = find(C==1); 

% combine values to a vector for reshaping
% in heritability variable 1st column A, 2nd column C, 3rd column E
W(Qind) = 0; W(Wind) = heritability(:,1); W(Eind) = 0; 
% reshape
heritMatrix = reshape(W,numNodes,numNodes); 
% make a full matrix
heritMatrix = heritMatrix+heritMatrix'; 

% make a curve plot
RichClubHuman(groupAdjlog,heritMatrix)

