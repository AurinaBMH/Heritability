%% for deterministic
load('GenCog_HCPMMP1ANDfslatlas20_FACT_connectome_data.mat')
load('GC_movement_params.mat')
type = standard;
length = standard_length;
numNodes = size(type{1},1);
Rmotion = zeros(numNodes);

threshold = 0.6;
selectConsist = false;
%% select only subjects that have connectomes
[subjects,ia] = intersect(subs, usedsubs);
m_rot = m_rot(ia)';
m_vol = m_vol(ia)';

motion = 'translation';
switch motion
    case 'translation'
        moparam = m_vol;
    case 'rotation'
        moparam = m_rot;
end

%% make 3D matrix from all connectomes
matrices = zeros(size(type,2), size(type{1},1),size(type{1},1));
avLength = zeros(size(type,2), size(type{1},1),size(type{1},1));
for subj=1:size(subjects,1)
    matrices(subj,:,:) = type{subj};
    avLength(subj,:,:) = length{subj};
end

% replace zeros with NaNs
avLength(avLength == 0) = NaN;
%% select the most consistent edges
adjs = logical(matrices);
nanMatrix = tril(nan(numNodes));
consist = triu(squeeze(mean(adjs,1)));
consist = consist+nanMatrix;
if selectConsist
    consistLog = consist>threshold;
    [row,col] = find(consistLog);
    ind = [row, col];
    for i=1:size(ind,1)

        weights = matrices(:,ind(i,1),ind(i,2));
        Rmotion(ind(i,1),ind(i,2)) = corr(weights,moparam);

    end
    consist = consist.*consistLog;

else
    for i=1:numNodes
        for j=i+1:numNodes

            weights = matrices(:,i,j);
            Rmotion(i,j) = corr(weights,moparam);

        end
    end

end

meanLength = squeeze(nanmean(avLength,1));
consist(consist == 0) = NaN;
Rmotion(Rmotion == 0) = NaN;
meanLength(isnan(Rmotion)) = NaN;

% calculate the probabilistic edge consistency
x = consist(~isnan(consist));
y = Rmotion(~isnan(Rmotion));
sz = 35;
c = meanLength(~isnan(meanLength));
%sz = lengthlinspace(1,100,200);
%scatter(x,y,sz)
figure; scatter(x, y, sz, c, 'filled'); title('Streamline count');
          lsline;
xlabel('Deterministic edge consistency (prop subj)'); ylabel('Motion effect on edge weight (r)')
[r,p] = corr(x,y);

% another way to plot the same data
figure; scatter(x,c, sz, y, 'filled');
xlabel('Deterministic edge consistency (prop subj)'); ylabel('Length (mm)');



%% for probabilistic
clear all;
load('GenCog_HCPMMP1ANDfslatlas20_iFOD2_connectome_data.mat')
load('GC_movement_params.mat')

type = standard_den;
numNodes = size(type{1},1);
Rmotion = zeros(numNodes);
%% select only subjects that have connectomes
[subjects,ia] = intersect(subs, usedsubs);
m_rot = m_rot(ia)';
m_vol = m_vol(ia)';

motion = 'rotation';
switch motion
    case 'translation'
        moparam = m_vol;
    case 'rotation'
        moparam = m_rot;
end

%% make 3D matrix from all connectomes
matrices = zeros(length(ia), size(type{1},1),size(type{1},1));
for subj=1:length(subjects)
    matrices(subj,:,:) = type{subj};
end
stdev = squeeze(std(matrices,1));
meanW = squeeze(mean(matrices,1));
coefVar = stdev./meanW;
consist = triu(-log(coefVar));
nanMatrix = tril(nan(numNodes));
consist = consist+nanMatrix;

for i=1:numNodes
    for j=i+1:numNodes

        weights = matrices(:,i,j);
        Rmotion(i,j) = corr(weights,moparam);

    end
end

Rmotion(Rmotion == 0) = NaN;

% calculate the probabilistic edge consistency
consist(~isfinite(consist)) = NaN;
x = consist(~isnan(consist));
y = Rmotion(~isnan(Rmotion));
figure; scatter(x, y, 'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5); xlabel('Probabilistic edge consistency -log(CoefVar))'); ylabel('Motion effect on edge weight (r)');
          title('Normalised streamline count');
[r,p] = corr(x, y);
