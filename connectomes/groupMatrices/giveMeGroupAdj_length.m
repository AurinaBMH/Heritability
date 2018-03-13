parcellation = 'HCPMMP1'; 
tract = 'FACT'; 
sift = 'SIFT2'; 
numBins = 20; 

C = load('%sANDfslatlas20_acpc_%s_%s_standard_structnets.mat', parcellation, tract, sift); 
L = load('%s_acpc_%s_%s_length_structnets.mat', parcellation, tract, sift); 

[numNodes, leftCortex, leftSubcortex, rightCortex, rightSubcortex] = labelNodes(parcellation); 

% make masks for within and between hemisheres
mask = zeros(NumNodes); 
maskWithin = zeros(NumNodes); 
maskWithin(min(LeftCortex):max(LeftSubcortex),min(LeftCortex):max(LeftSubcortex)) = 1; 
maskWithin(min(RightCortex):max(RightSubcortex),min(RightCortex):max(RightSubcortex)) = 1;
maskBetween = ~maskWithin; 

% make 3D matrix and calculate mean distances
length = zeros(size(C.ADJS{1},1),size(C.ADJS{1},1), size(C.ADJS,2)); 
adj = zeros(size(C.ADJS{1},1),size(C.ADJS{1},1), size(C.ADJS,2)); 
consist = zeros(size(C.ADJS{1},1),size(C.ADJS{1},1), size(C.ADJS,2)); 

for s = 1:size(C.ADJS,2)
    length(:,:,s) = L.ADJS{s}; 
    adj(:,:,s) = C.ADJS{s}; 
    consist(:,:,s) = logical(C.ADJS{s}); 
end

mLength = mean(length,3); 
mconsist = mean(consist,3); 


lengthWithin = mLength.*maskWithin; 
lengthWithin(lengthWithin==0) = NaN; 
thresholdsWithin = quantile(lengthWithin(:),linspace(0,1,numBins+1));
binWithin = discretize(lengthWithin,thresholdsWithin, 'IncludedEdge','right'); 

lengthBetween = mLength.*maskBetween; 
lengthBetween(lengthBetween==0) = NaN; 
thresholdsBetween = quantile(lengthBetween(:),linspace(0,1,numBins+1));
binBetween = discretize(lengthBetween,thresholdsBetween, 'IncludedEdge','right'); 





