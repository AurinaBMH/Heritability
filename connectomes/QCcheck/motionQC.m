
% the effect of motion
clear all; close all; 

parcellation = 'custom500'; % 'HCPMMP1' , 'custom200';
tract = 'FACT';
sift = 'SIFT2';
weight = 'FA';% 'strandard'
doPlot = true;
threshold = 0.5;
everySecond = false;


cd ('data/connectomes')
Conn = load(sprintf('%sANDfslatlas20_acpc_%s_%s_%s_structnets.mat', parcellation, tract, sift, weight));
Length = load(sprintf('%sANDfslatlas20_acpc_%s_%s_length_structnets.mat', parcellation, tract, sift));
M = load('motion_parameters.mat');

% calculate the average on b0 volumes
load('Volume_index.mat')
b0 = [1 17 33 49 65 81 96 113 129 145 161 177 192 209 225 241 257 273];
[b0IND, a] = find(volume_ind==b0);

coordinates = Conn.COG;
connectomes = Conn.ADJS;
distances = Length.ADJS;
motion = M.motion;

numNodes = size(coordinates{1},1);
numSubj = size(coordinates,2);
numEdges = (numNodes*numNodes-numNodes)/2;

% calculate the average of motion
avMotion = zeros(numSubj,2);
adj = zeros(numNodes, numNodes, numSubj);
avMotionb0 = zeros(numSubj,1);
secondIND = 1:2:576;
secondINDb0 = 1:2:36;
deg = zeros(numSubj,numNodes);
skew = zeros(numSubj,1);
density = zeros(numSubj,1);
strength = zeros(numSubj,1);
binAdj = zeros(numNodes, numNodes, numSubj);
maskAdj = zeros(numNodes, numNodes, numSubj);
edgeLength = zeros(numNodes, numNodes, numSubj);
edgesAdj = zeros(numSubj,numEdges);

for s=1:numSubj
    
    b0motion = motion{s}(b0IND,1);
    
    if everySecond
        allmotion = motion{s}(secondIND,:);
        b0motion = b0motion(secondINDb0);
    else
        allmotion = motion{s};
    end
    avMotionb0(s) = mean(b0motion);
    avMotion(s,:) = mean(allmotion,1);
    adj(:,:,s) = connectomes{s};
    % calculate skewness for each subject to use it as "hubness"
    deg(s,:) = degrees_und(connectomes{s});
    density(s) = density_und(connectomes{s});
    skew(s) = skewness(deg(s,:));
    strength(s) = sum(nonzeros(connectomes{s}));
    binAdj(:,:,s) = logical(connectomes{s});
    maskAdj(:,:,s) = maskuHalf(adj(:,:,s));
    mask = ~isnan(maskAdj(:,:,s));
    adjSubj = maskAdj(:,:,s);
    edgesAdj(s,:) = adjSubj(mask==1);
    edgeLength(:,:,s) = distances{s};
    
end

b0motion = avMotion(b0IND,1);
%figure; histogram(avMotion(:,1), 50); xlabel('motion - all volumes');
%figure; histogram(avMotionb0, 50); xlabel('motion - b0 volumes');

% HOW CORRELATED MOTION ON B0 AND ALL VOLUMES;
[r] = corr(avMotion(:,1), avMotionb0, 'type', 'Spearman');
figure; scatter(avMotion(:,1), avMotionb0);
xlabel('motion all volumes'); ylabel('motion b0 volumes')

% corelate skewness with motion
rskewb0 = corr(skew, avMotionb0, 'type', 'Spearman');
rskew = corr(skew, avMotion(:,1), 'type', 'Spearman');

%figure; scatter(avMotionb0, skew); xlabel('motion - b0 volumes');  ylabel('Skewness');
%figure; scatter(avMotion(:,1), skew); xlabel('motion - all volumes');  ylabel('Skewness');

% correlate density with motion
rdensb0 = corr(density, avMotionb0, 'type', 'Spearman');
rdens = corr(density, avMotion(:,1), 'type', 'Spearman');

%figure; scatter(avMotionb0, density); xlabel('motion - b0 volumes');  ylabel('Density');
%figure; scatter(avMotion(:,1), density); xlabel('motion - all volumes');  ylabel('Density');

% correlate total strength of the connectome with motion
rstrengthb0 = corr(strength, avMotionb0, 'type', 'Spearman');
rstrength = corr(strength, avMotion(:,1), 'type', 'Spearman');

%figure; scatter(avMotionb0, strength); xlabel('motion - b0 volumes');  ylabel('Strength');
%figure; scatter(avMotion(:,1), strength); xlabel('motion - all volumes');  ylabel('Strength');

% reproduce baum


if strcmp(tract, 'iFOD2')
    [~, consistency] = giveMeGroupAdj_variance(connectomes);
    consistency(isnan(consistency)) = 0; 
    consistency = maskuHalf(consistency);
    consistency(isnan(consistency)) = [];
    consistency = -log(consistency); 
    indKeep = find(isfinite(consistency)); 
    %consistency = consistency(indKeep); 
    %consistency(~isfinite(consistency))=0; 
 
elseif strcmp(tract, 'FACT')
    consistency = mean(binAdj,3);
    consistency = maskuHalf(consistency);
    consistency(isnan(consistency)) = [];
    indKeep = find(consistency~=0);
end

% % for each edge correlate consistency with motion
% rcons = zeros(length(consistency),1);
% rconsR = zeros(length(consistency),1);
%
% for edge = 1:length(consistency)
%     rcons(edge) = corr(consistency, avMotion);
%     rconsR(edge) = corr(consistencyRoberts, avMotion);
% end
%
% figure; histogram(rcons, 100);
% figure; histogram(rconsR, 100);

% for each edge correlate weight with motion
rweight = zeros(length(consistency),1);
pweight = zeros(length(consistency),1);
for edge = 1:numEdges
    
    [rweight(edge), pweight(edge)] = corr(edgesAdj(:,edge), avMotion(:,1), 'type', 'Spearman');
    
end

% find proportion of significant correlations
thresh = 0.05; 
P = sum(pweight(indKeep)<thresh)/length(indKeep);
fprintf('%d edges are signifficantly affected by motion at p=%d\n', P, thresh)
rweight(isnan(rweight)) = 0;
figure; histogram(rweight(indKeep), 100); xlabel('correlation between weight and motion'); ylabel('number of edges');
%figure; histogram(pweight, 100); xlabel('correlation between weight and motion'); ylabel('number of edges');

% do with length
% calulate length of each edge
avLength = sum(edgeLength,3) ./ sum(edgeLength~=0,3);
avLength(isnan(avLength)) = 0;
avLength = maskuHalf(avLength);
avLength = avLength(~isnan(avLength));

nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
% bin length into a hundred bins
numBins = 100; 
thresholdsLength = quantile(avLength(indKeep),linspace(0,1,numBins+1));
binLength = discretize(avLength(indKeep),thresholdsLength, 'IncludedEdge','right'); 

%binLength = discretize(avLength,100);
cLength = zeros(length(binLength),3);
for i=1:length(binLength)
    
    cLength(i,:) = nice_cmap(binLength(i),:);
end
figure; sz = 30;
indKeep2 = find(consistency~=0); 
scatter(consistency(indKeep)', rweight(indKeep), sz,cLength, 'filled'); 
xlabel('Consistency'); ylabel('Motion effect r'); title ('Length in color');

% bin r into a hundred bins
thresholdsR = quantile(rweight(indKeep),linspace(0,1,numBins+1));
binR = discretize(rweight(indKeep),thresholdsR, 'IncludedEdge','right'); 

%binR = discretize(rweight,100);
cR = zeros(length(binR),3);
for i=1:length(binLength)
    cR(i,:) = nice_cmap(binR(i),:);
end

figure; sz = 30;
scatter(consistency(indKeep)', avLength(indKeep), sz,cR, 'filled'); 
xlabel('Consistency'); ylabel('Length'); title ('Motion effect (r) in color');