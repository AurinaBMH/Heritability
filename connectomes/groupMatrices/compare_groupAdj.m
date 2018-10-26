% make average matrices
parcellation = 'HCPMMP1'; %'HCPMMP1' , 'custom200'
tract = 'iFOD2';
sift = 'SIFT2';
weight = 'standard'; 
threshold = 0.3;

Conn = load(sprintf('%sANDfslatlas20_acpc_%s_%s_%s_structnets.mat', parcellation, tract, sift, weight));
Length = load(sprintf('%sANDfslatlas20_acpc_%s_%s_length_structnets.mat', parcellation, tract, sift));

coordinates = Conn.COG;
connectomes = Conn.ADJS;
distances = Length.ADJS;

% calculate average distances between ROIs as an average of all subjects
numNodes = size(coordinates{1},1);
numSubj = size(coordinates,2);

% make vectors for hemispheres for different parcellations based on the
% number of nodes (first half is always left, second half is always right)
hemiid = zeros(numNodes,1);
hemiid(1:numNodes/2) = 1;
hemiid(numNodes/2+1:numNodes) = 2;

dist = zeros(numNodes, numNodes, numSubj);
adjMatr = zeros(numNodes, numNodes, numSubj);

for s=1:numSubj
    
    dist(:,:,s) = pdist2(coordinates{s}, coordinates{s});
    adjMatr(:,:,s) = connectomes{s};
end
% take average of distance
avDist = mean(dist,3);

[groupAdj_variance, consist_var] = giveMeGroupAdj_variance(connectomes, threshold);
[groupAdj_consistency, groupDist, consist_cons] = giveMeGroupAdj_consistency(connectomes, distances,threshold);
% distance bin consistency (sent by Bratislav Misic) - this gives binary
% group matrix; can use that as mask and take the average of nonzero values
% across subjects for those edges
groupAdj_maskLength = fcn_group_average(adjMatr,avDist,hemiid);
% replace zeros with NaNs; 
adjMatr(adjMatr==0) = NaN; 
% take average across subjectsclose all
meanAdj = nanmean(adjMatr,3); 
groupAdj_length = meanAdj.*groupAdj_maskLength; 
groupAdj_length(isnan(groupAdj_length)) = 0; 

% compare resulting group matrices by plotting and calculating correlations
figure; scatter(groupAdj_consistency(:), groupAdj_variance(:)); xlabel('consistency'); ylabel('variance'); 
r1 = corr(groupAdj_consistency(:), groupAdj_variance(:), 'type', 'spearman'); 

figure; scatter(groupAdj_consistency(:), groupAdj_length(:)); xlabel('consistency'); ylabel('length'); 
r2 = corr(groupAdj_consistency(:), groupAdj_length(:), 'type', 'spearman'); 

figure; scatter(groupAdj_variance(:), groupAdj_length(:)); xlabel('variance'); ylabel('length'); 
r3 = corr(groupAdj_variance(:), groupAdj_length(:), 'type', 'spearman');
