
function [rall, pall] = connectomeCorrelation(connectomes, doPlot)


% calculate subject-subject correlations
num_datasets = length(connectomes);
r = zeros(num_datasets);
p = zeros(num_datasets);

for s1 = 1:num_datasets
    for s2 = s1+1:num_datasets
        % take the lower half of the matrix
        half1 = maskuHalf(connectomes{s1});
        half2 = maskuHalf(connectomes{s2});
        % select existing values into a vector
        weights1 = half1(~isnan(half1));
        weights2 = half2(~isnan(half2));
        % correlate each pair of subjects
        [r(s1,s2),p(s1,s2)] = corr(weights1, weights2, 'type', 'Spearman');

    end
end

rall = r+r';
pall = p+p';

if doPlot
    figure; imagesc(rall);
    set(gcf,'color','w');
    nice_cmap = [make_cmap('steelblue',50,30,0);flipud(make_cmap('orangered',50,30,0))];
    colormap(nice_cmap)
    caxis([0 1])
end
end
% check matrices for subjects that have correlation <0.6 to others.
% subject 299 shows consistently low correlation with all other subjects -
% low number of inter-hemispheric connections

% compare scatterplots between regular and 299 subjects
