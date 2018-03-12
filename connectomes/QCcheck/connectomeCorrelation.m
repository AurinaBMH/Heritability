cd ('data/connectomes')
load('HCPMMP1ANDfslatlas20_acpc_FACT_SIFT2_standard_structnets.mat')

% make an average histogram for all subjects; 
edges = 1:2:120; 
deg = zeros(length(ADJS), size(ADJS{1},1));
numedges = length(edges);
num_datasets = length(ADJS);
counts = zeros(numedges-1, num_datasets);



for dataset_idx = 1:num_datasets
    deg(dataset_idx,:) = degrees_und(ADJS{dataset_idx});
    counts(:, dataset_idx) = histcounts(deg(dataset_idx,:), edges);
end

avg_counts = mean(counts,2);
figure; 
bar(edges(1:end-1), avg_counts); title(sprintf('Average of histograms for %d subjects', num_datasets))

% calculate subject-subject correlations
r = zeros(num_datasets); 
p = zeros(num_datasets); 

for s1 = 1:num_datasets
    for s2 = s1+1:num_datasets
        % take the lower half of the matrix
        half1 = maskuHalf(ADJS{s1}); 
        half2 = maskuHalf(ADJS{s2}); 
        % select existing values into a vector
        weights1 = half1(~isnan(half1)); 
        weights2 = half2(~isnan(half2));
        % correlate each pair of subjects
        [r(s1,s2),p(s1,s2)] = corr(weights1, weights2, 'type', 'Spearman'); 
        
    end
end

rall = r+r'; 
meanR = mean(rall,1); 
% check matrices for subjects that have correlation <0.6 to others. 
% subject 299 shows consistently low correlation with all other subjects -
% low number of inter-hemispheric connections

% compare scatterplots between regular and 299 subjects
figure; scatter(ADJS{298}(:), ADJS{300}(:)); hold on; 
scatter(ADJS{299}(:), ADJS{301}(:));
legend({'298 and 300', '299 and 301'})

