% make average matrices
parcellation = 'HCPMMP1'; %'HCPMMP1' , 'custom200'
tract = 'iFOD2';
sift = 'SIFT2';
thresholds = [0.3 0.15 0.05];
weight = 'standard';
calculateRC = false;
onlyCortex = true; 

% plotting options
scatterColor = [.98 .85 .37];
histoColor = [.28 .75 .57];
edgeColor = [0.45 0.45 0.45];

Conn = load(sprintf('%sANDfslatlas20_acpc_%s_%s_%s_structnets.mat', parcellation, tract, sift, weight));
Length = load(sprintf('%sANDfslatlas20_acpc_%s_%s_length_structnets.mat', parcellation, tract, sift));

coordinates = Conn.COG;
connectomes = Conn.ADJS;
distances = Length.ADJS;

% calculate average distances between ROIs as an average of all subjects
if onlyCortex
    numNodes = size(coordinates{1},1)-20;
else
    
    numNodes = size(coordinates{1},1);
end
numSubj = size(coordinates,2);

% make vectors for hemispheres for different parcellations based on the
% number of nodes (first half is always left, second half is always right)
hemiid = zeros(numNodes,1);
hemiid(1:numNodes/2) = 1;
hemiid(numNodes/2+1:numNodes) = 2;

dist = zeros(numNodes, numNodes, numSubj);
adjMatr = zeros(numNodes, numNodes, numSubj);

for s=1:numSubj
    
    if onlyCortex
    dist(:,:,s) = pdist2(coordinates{s}([1:numNodes/2, numNodes/2+11:numNodes+10], :), coordinates{s}([1:numNodes/2, numNodes/2+11:numNodes+10],:));
    adjMatr(:,:,s) = connectomes{s}([1:numNodes/2, numNodes/2+11:numNodes+10], [1:numNodes/2, numNodes/2+11:numNodes+10]);
    connectomes{s} = connectomes{s}([1:numNodes/2, numNodes/2+11:numNodes+10], [1:numNodes/2, numNodes/2+11:numNodes+10]);
    else
    dist(:,:,s) = pdist2(coordinates{s}, coordinates{s});
    adjMatr(:,:,s) = connectomes{s};
    end
end
% take average of distance
avDist = mean(dist,3);
k=1;
%figure;
for thr=thresholds
    [groupAdj_variance] = giveMeGroupAdj_variance(connectomes, thr);
    
    % plot the degree distribution for the group matrix
    [deg] = degrees_und(groupAdj_variance); strength = strengths_und(groupAdj_variance);
    
    figure;
    %subplot(3,2,k);
    histogram(deg, 50, 'EdgeColor', edgeColor, 'FaceColor', histoColor);
    title (sprintf('degree distribution %s%', thr)); xlabel('Degree')
    hold on;
    
    % plot the strength distribution of the group matrix
    %subplot(3,2,k+3);
    figure;
    histogram(strength, 50, 'EdgeColor', edgeColor, 'FaceColor', histoColor);
    title (sprintf('Strength distribution %s%', thr)); xlabel('Strength')
    hold on;
    
    if calculateRC
        % calculate topological RC
        numIter = 50;
        numRepeats = 100;
        WhatTypeNetwork = 'bu';
        whatNullModel = 'randmio_und'; %'randmio_und'; %'strength'; %
        
        if strcmp(WhatTypeNetwork, 'bu')
            groupAdj_variance = logical(groupAdj_variance);
        end
        
        if strcmp(WhatTypeNetwork, 'wuStrength') || strcmp(WhatTypeNetwork, 'wuStrengthBins')
            kmax = 100;
        else
            kmax = max(sum(logical(groupAdj_variance)));
        end
        strength = sum(groupAdj_variance);
        deg = degrees_und(groupAdj_variance);
        pThreshold = 0.05;
        %[~, E] = BF_PlotQuantiles(strength,strength,100);
        [PhiNorm,PhiTrue,PhiRand] = RichClubPhiNorm(groupAdj_variance,kmax, numIter,numRepeats,WhatTypeNetwork,whatNullModel); %, doBins);
        figure;
        % Compute p-values
        pValues = zeros(kmax,1);
        for i = 1:kmax
            pValues(i) = mean(PhiTrue(i) <= PhiRand(:,i));
        end
        % Significant phi
        isSig = (pValues <= pThreshold);
        PhiNormMean = zeros(size(PhiTrue));
        for i = 1:length(PhiTrue)
            PhiNormMean(i) = PhiTrue(i)/mean(PhiRand(:,i));
        end
        % plot the graphs
        if strcmp(WhatTypeNetwork, 'wuStrength') || strcmp(WhatTypeNetwork, 'wuStrengthBins')
            
            subplot(2,1,1); plot(PhiNormMean, '-','Color','r','LineWidth',2);  xlim([0 100]); ylim([min(PhiNormMean)-0.1 max(PhiNormMean)+0.05]); legend('Phi normalised', 'Location','NorthWest'); ...
                hold on;
            plot(find(isSig),PhiNormMean(isSig),'o','Color','r','LineWidth',3); xlim([0 100]); ylim([min(PhiNormMean)-0.1 max(PhiNormMean)+0.05]);
            title (sprintf('Normalised rich club \n%s', type)); xlabel('Strength'); ylabel('Phi (norm)'); set(gca,'FontSize',12,'fontWeight','bold');
            subplot(2,1,2); histogram(strength, 50, 'EdgeColor', edgeColor, 'FaceColor', histoColor);  title ('Strength distribution');  xlabel('Strength bins'); set(gca,'FontSize',12,'fontWeight','bold');
            %     else
            %histogram(strength, 50);
            %         subplot(2,1,1); plot(PhiNormMean, '-','Color','r','LineWidth',2);  xlim([0 100]);
            %         ylim([min(PhiNormMean)-0.1 max(PhiNormMean)+0.05]); legend('Phi normalised', 'Location','NorthWest'); ...
            %         hold on;
            %         plot(find(isSig),PhiNormMean(isSig),'o','Color','r','LineWidth',3); xlim([min(E) max(E)]); ylim([min(PhiNormMean)-0.1 max(PhiNormMean)+0.05]);
            %         title (sprintf('Normalised rich club \n%s', type)); xlabel('Strength'); ylabel('Phi (norm)'); set(gca,'FontSize',12,'fontWeight','bold');
            %subplot(2,1,2); histogram(deg, 50); title ('Degree distribution'); xlim([min(deg) max(deg)+2]); xlabel('Degree'); set(gca,'FontSize',12,'fontWeight','bold');
            
        else
            
            subplot(2,1,1); plot(PhiNormMean, '-','Color','r','LineWidth',2);  xlim([min(deg) max(deg)+2]); ylim([0.8 max(PhiNormMean)+0.05]); legend('Phi normalised', 'Location','NorthWest'); ...
                hold on;
            plot(find(isSig),PhiNormMean(isSig),'o','Color','r','LineWidth',3); xlim([min(deg) max(deg)+2]); ylim([0.8 max(PhiNormMean)+0.05]);
            title (sprintf('Normalised rich club')); xlabel('Degree'); ylabel('Phi (norm)'); set(gca,'FontSize',12,'fontWeight','bold');
            subplot(2,1,2); histogram(deg, 50, 'EdgeColor', edgeColor, 'FaceColor', histoColor); title ('Degree distribution'); xlim([min(deg) max(deg)+2]); xlabel('Degree'); set(gca,'FontSize',12,'fontWeight','bold');
        end
    end
    k=k+1;
end
