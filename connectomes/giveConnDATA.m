
function [A, matrices, coordinates, avWeight, SUBS] = giveConnDATA(whatDWI,parc,tract,weight,brainPart,nRem)

if nargin<5
    weight = 'standard';
    brainPart = 'Lcortex';
end


if strcmp(parc, 'aparcaseg')
    parcellation = 'aparc+aseg';
    LC = 1:34;
    RC = 42:75;
    nROIs = 82;
elseif strcmp(parc, 'HCP')
    parcellation = 'HCPMMP1ANDfslatlas20';
    LC = 1:180;
    RC = 191:370;
    nROIs = 360;
elseif strcmp(parc, 'cust100')
    parcellation = 'custom200ANDfslatlas20';
    LC = 1:100;
    RC = 111:210;
    nROIs = 220;
elseif strcmp(parc, 'cust250')
    parcellation = 'custom500ANDfslatlas20';
    LC = 1:250;
    RC = 261:510;
    nROIs = 530;
end

switch whatDWI
    case 'HCP'
        load(sprintf('%s_acpc_%s_NOSIFT_%s_structnets.mat', parcellation, tract, weight));
        % remove subject 299
        ADJS{299} = []; ADJS = ADJS(~cellfun('isempty',ADJS));
        COG{299} = []; COG = COG(~cellfun('isempty',COG));
        SUBS(299) = [];
    case 'GenCog'
        load('GenCog_combinations_iFOD2_FACT_NOS_onlyTEST.mat')
        ADJS = connectomes.(tract).newACT.dynamic.noSIFT.NOS.(parc);
        ADJS = ADJS(~cellfun('isempty',ADJS));
        COG = connectomes.COG.(parc);
end

cort = [LC,RC];
if strcmp(brainPart, 'LRcortex')
    keepNodes = cort;
elseif strcmp(brainPart, 'Lcortex')
    keepNodes = LC;
elseif strcmp(brainPart, 'wholeBrain')
    keepNodes = 1:size(ADJS{1},1);
end

numNodes = length(keepNodes);
numSubjects = length(ADJS);
matrices = cell(1,size(ADJS,2));
coordinates = cell(1,size(ADJS,2));
A = zeros(numNodes,numNodes,numSubjects);

for m=1:length(ADJS)
    % remove 25% of weakest links
    ADJS{m}(ADJS{m}<nRem) = 0; % remove less than 10 streamlines
    if ~strcmp(brainPart, 'wholeBrain')
        matrices{m} = ADJS{m}(keepNodes, keepNodes);
        A(:,:,m) = matrices{m};
        switch whatDWI
            case 'HCP'
                coordinates{m} = COG{m}(keepNodes,:);
            case 'GenCog'
                coordinates{m} = COG(keepNodes,:); % assigne the same coords for all subjects
        end
        
    else
        matrices{m} = ADJS{m};
        A(:,:,m) = matrices{m};
        switch whatDWI
            case 'HCP'
                coordinates{m} = COG{m};
            case 'GenCog'
                coordinates{m} = COG; % assigne the same coords for all subjects
        end
    end
    
end

nanA = A;
nanA(nanA==0) = NaN;
avWeight = nanmean(nanA,3);
avWeight(isnan(avWeight)) = 0;
end