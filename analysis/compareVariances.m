ConnMask = load('HCPMMP1ANDfslatlas20_acpc_iFOD2_SIFT2_standard_structnets.mat');
ConnLength = load('HCPMMP1ANDfslatlas20_acpc_iFOD2_SIFT2_length_structnets.mat');  

load('twinEdges_HCPMMP1_iFOD2_FATEST15.mat')
load('twinCovariatesDWI.mat')

varianceA = zeros(size(Output_MZ,3),1);
varianceB = zeros(size(Output_MZ,3),1);
varianceT1 = zeros(size(Output_MZ,3),1);
varianceT2 = zeros(size(Output_MZ,3),1);
for i=1:size(Output_MZ,3)
    
A1 = Output_DZ(:,:,i);
A1(A1==0) = NaN; 
B1 = Output_MZ(:,:,i);
B1(B1==0) = NaN;

A2 = A1(:);
B2 = B1(:);
% variance by zigosity
varianceA(i) = nanvar(A2);
varianceB(i) = nanvar(B2);

% variance by twin1/2
T1 = vertcat(A1(:,1), B1(:,1));
T2 = vertcat(A1(:,2), B1(:,2));

varianceT1(i) = nanvar(T1);
varianceT2(i) = nanvar(T2);
end

numNodes = size(ConnMask.ADJS{1},1);
twins = vertcat(MZ_ID(:), DZ_ID(:));
[~, twinIND] = intersect(ConnMask.SUBS, twins);

coordinatesMask = ConnMask.COG(twinIND);
connectomesMask = ConnMask.ADJS(twinIND);

coordinatesLength = ConnLength.COG(twinIND);
connectomesLength = ConnLength.ADJS(twinIND);

[groupAdj, consist_var] = giveMeGroupAdj_variance(connectomesMask, 0.15);
groupAdjlog = logical(groupAdj); 

groupAdjLength = ConnLength.ADJS{2};
groupAdjLength = maskuHalf(groupAdjLength);
groupAdjLength(groupAdjlog==0) = NaN;

edgeLabelLength = groupAdjLength(~isnan(groupAdjLength));
figure; histogram(edgeLabelLength, 50);

deg = degrees_und(groupAdjlog);

colorHub = false; 

for k=80
    isHub = deg>k;
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
    %get indeces for existing values in group matrix (just on half so they can
    %be used to ut heritability values back to the matrix)
    cols = zeros(size(edgeLabel,1),3);
    for j=1:size(edgeLabel,1)
        if colorHub
            if edgeLabel(j)==3
                cols(j,:) = [.89 0 .06];
            elseif edgeLabel(j)==2
                cols(j,:) = [.97 .81 .16];
            elseif edgeLabel(j)==1
                cols(j,:) = [0 .47 .75]; %
            end
        else
            if edgeLabelLength(j)<80
                cols(j,:) = [0 .47 .75];
            elseif edgeLabelLength(j)>=80 && edgeLabelLength(j)<150
                cols(j,:) = [.97 .81 .16];
            elseif edgeLabelLength(j)>=150
                cols(j,:) = [.89 0 .06];
            end
        end
        
    end
    figure;
        subplot(1,2,1); scatter(varianceA/100000000, varianceB/100000000, 20, cols, 'filled'); ...
            xlabel('DZ'); ylabel('MZ'); title(sprintf('Hub threshold is %d', k));
        subplot(1,2,2); scatter(varianceT1/100000000, varianceT2/100000000, 20, cols, 'filled'); ...
            xlabel('twin1'); ylabel('twin2'); title(sprintf('Hub threshold is %d', k));
end