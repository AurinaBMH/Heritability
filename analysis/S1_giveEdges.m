function S1_giveEdges(whatDATA,parcellation,tractography,brainPart,weight1,strRem1,densThreshold, groupConn,cvMeasure,consThr,weight2,strRem2)

[A, matrices, coordinates, avWeight, SUBjects] = giveConnDATA(whatDATA,parcellation,tractography,weight1,brainPart,strRem1); 
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

giveRC=false; 


[Gr] = giveMeRichClub(matrices, coordinates, groupConn ,densThreshold, giveRC, cvMeasure, consThr);
groupAdjlog = logical(Gr); 
dens = density_und(Gr); 

%-----------------------------------------------------------------
% Select edges using FA as weight
%-----------------------------------------------------------------

[A, matrices, coordinates, avWeight, SUBjects] = giveConnDATA(whatDATA,parcellation,tractography,weight2,brainPart,strRem2); 
coordinatesMask = coordinates(subIND);
connectomesMask = matrices(subIND);

% get indeces for non-zero elements

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

pathname = fileparts('data/output/');
%use that when you save
matfileEdges = fullfile(pathname, sprintf('twinEdges_%s_%s_%s_%s%d.mat', parcellation, tractography, weight2, cvMeasure, round(densThreshold*100)));
save(matfileEdges, 'Output_MZ', 'Output_DZ', 'groupAdjlog');

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
end