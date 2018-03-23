% make a txt file for an edge
clear all; close all;
parcellation = 'custom200'; % 'HCPMMP1' , 'custom200';
tract = 'FACT';
tract2 = 'iFOD2';
sift = 'SIFT2';
groupConn = 'length';
weight = 'standard'; % 'FA', 'standard'

cd ('data/general')

load('twinCovariatesDWI.mat')
cd ..
cd ('connectomes')

ConnMask = load(sprintf('%sANDfslatlas20_acpc_%s_%s_%s_structnets.mat', parcellation, tract, sift, weight));

twins = vertcat(MZ_ID(:), DZ_ID(:));
[~, twinIND] = intersect(ConnMask.SUBS, twins);

coordinatesMask = ConnMask.COG(twinIND);
connectomesMask = ConnMask.ADJS(twinIND);


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
    [groupAdj, consist_var] = giveMeGroupAdj_variance(connectomesMask);
elseif strcmp(groupConn, 'consistency')
    
    [groupAdj, groupDist, consist_cons] = giveMeGroupAdj_consistency(connectomesMask, distances,threshold);
end

groupAdjlog = logical(groupAdj);

ConnHerit = load(sprintf('%sANDfslatlas20_acpc_%s_%s_%s_structnets.mat', parcellation, tract2, sift, weight));

[~, twinIND] = intersect(ConnHerit.SUBS, twins);

coordinatesHerit = ConnHerit.COG(twinIND);
connectomesHerit = ConnHerit.ADJS(twinIND);
SUBS = ConnHerit.SUBS(twinIND);


edgeMatr = zeros(sum(groupAdjlog(:))/2, numSubj);
for s=1:numSubj
    m = maskuHalf(connectomesHerit{s}.*groupAdjlog);
    m(groupAdjlog==0) = NaN;
    edgeMatr(:,s) = m(~isnan(m));
end

numEdges = size(edgeMatr,1);
Output_MZ = zeros(size(MZ_ID,1),2, numEdges);

for i=1:size(MZ_ID,1)
    ind1 = find(SUBS==MZ_ID(i,1));
    ind2 = find(SUBS==MZ_ID(i,2));
    for edg = 1:numEdges
        Output_MZ(i,1,edg) = edgeMatr(edg,ind1);
        Output_MZ(i,2,edg) = edgeMatr(edg,ind2);
    end
end

Output_DZ = zeros(size(DZ_ID,1),2, numEdges);
for i=1:size(DZ_ID,1)
    ind1 = find(SUBS==DZ_ID(i,1));
    ind2 = find(SUBS==DZ_ID(i,2));
    for edg = 1:numEdges
        Output_DZ(i,1,edg) = edgeMatr(edg,ind1);
        Output_DZ(i,2,edg) = edgeMatr(edg,ind2);
    end
end

cd ..
cd ('output')
fileName = sprintf('twinEdges_%s_%s_%s.mat', parcellation, tract2, weight);
save(fileName, 'Output_MZ', 'Output_DZ');