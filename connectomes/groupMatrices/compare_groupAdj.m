% make average matrices
parcellation = 'HCPMMP1'; 
tract = 'FACT'; 
sift = 'SIFT2'; 
threshold = 0.5; 

C = load(sprintf('%sANDfslatlas20_acpc_%s_%s_standard_structnets.mat', parcellation, tract, sift)); 
L = load(sprintf('%sANDfslatlas20_acpc_%s_%s_length_structnets.mat', parcellation, tract, sift)); 

connectomes = C.ADJS; 
distances = L.ADJS; 

[groupAdj_variance, consist_var] = giveMeGroupAdj_variance(connectomes); 
[groupAdj_consistency, groupDist, consist_cons] = giveMeGroupAdj_consistency(connectomes, distances,threshold); 
% distance bin consistency
