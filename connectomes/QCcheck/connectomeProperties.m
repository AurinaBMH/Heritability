cd ('data/connectomes')
%% Define options
parcellation = 'custom200'; % 'HCPMMP1' , 'custom200'; 
tract = 'FACT';
sift = 'SIFT2';
doPlot = true; 
threshold = 0.5; 

%-----------------------------------------------------------------
% Arrange data 
%-----------------------------------------------------------------
% plot the histogram for the average matrix
Conn = load(sprintf('%sANDfslatlas20_acpc_%s_%s_standard_structnets.mat', parcellation, tract, sift));
Length = load(sprintf('%sANDfslatlas20_acpc_%s_%s_length_structnets.mat', parcellation, tract, sift));

coordinates = Conn.COG;
connectomes = Conn.ADJS;
distances = Length.ADJS;

%-----------------------------------------------------------------
% Plot the average of histograms
%-----------------------------------------------------------------
avg_counts = averageHistogram(connectomes, doPlot); 

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

%-----------------------------------------------------------------
% Make three versions of group matrices
%-----------------------------------------------------------------
% plot the histogram for the average matrix
% create group matrices
[groupAdj_variance, consist_var] = giveMeGroupAdj_variance(connectomes);
[groupAdj_consistency, groupDist, consist_cons] = giveMeGroupAdj_consistency(connectomes, distances,threshold);groupAdj_maskLength = fcn_group_average(adjMatr,avDist,hemiid);
% replace zeros with NaNs; 
groupAdj_maskLength = fcn_group_average(adjMatr,avDist,hemiid);
% replace zeros with NaNs; 
adjMatr(adjMatr==0) = NaN; 
% take average across subjectsclose all
meanAdj = nanmean(adjMatr,3); 
groupAdj_length = meanAdj.*groupAdj_maskLength; 
groupAdj_length(isnan(groupAdj_length)) = 0; 

%-----------------------------------------------------------------
% Plot degree distribution for each
%-----------------------------------------------------------------
deg_var = degrees_und(groupAdj_variance); 
deg_cons = degrees_und(groupAdj_consistency); 
deg_length = degrees_und(groupAdj_length); 

figure; 
subplot(3,1,1); histogram(deg_var, 30, 'FaceColor', [1 .43 .29], 'EdgeColor', [.45 .45 .45]); title('Degree distribution variance-based group matrix'); xlabel('Degree, k'); 
subplot(3,1,2); histogram(deg_cons, 30, 'FaceColor', [1 .43 .29], 'EdgeColor', [.45 .45 .45]); title('Degree distribution consistency-based group matrix'); xlabel('Degree, k'); 
subplot(3,1,3); histogram(deg_length, 30, 'FaceColor', [1 .43 .29], 'EdgeColor', [.45 .45 .45]); title('Degree distribution distance-based group matrix'); xlabel('Degree, k'); 

%-----------------------------------------------------------------
% Calculate RC curves for group matrices
%-----------------------------------------------------------------
WhatTypeNetwork = 'wu'; % 'wu' - weighted undirected; 'bu' - binary undirected; 
whatNullModel = 'randmio_und'; % 'randmio_und' - randomise topology; 'shuffleWeights' - keep topology, randomise weights. 

RCcurves(groupAdj_length, WhatTypeNetwork,whatNullModel)
title(sprintf('Group connectome - length %s - %s - %s %s %s',parcellation, tract, sift, WhatTypeNetwork, whatNullModel)); 
RCcurves(groupAdj_variance, WhatTypeNetwork,whatNullModel)
title(sprintf('Group connectome - variance %s - %s - %s %s %s',parcellation, tract, sift, WhatTypeNetwork, whatNullModel)); 
RCcurves(groupAdj_consistency, WhatTypeNetwork,whatNullModel)
title(sprintf('Group connectome - consistency %s - %s - %s %s %s',parcellation, tract, sift, WhatTypeNetwork, whatNullModel)); 

%-----------------------------------------------------------------
% Calculate RC curves for individuals
%-----------------------------------------------------------------


%-----------------------------------------------------------------
% Calculate correlations between individual connectomes
%-----------------------------------------------------------------
%[rall, pall] = connectomeCorrelation(parcellation, tract, sift, doPlot); 


