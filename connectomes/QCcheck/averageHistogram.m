function avg_counts = averageHistogram(connectomes, doPlot)

% make an average histogram for all subjects;
edges = 1:2:120;
deg = zeros(length(connectomes), size(connectomes{1},1));
numedges = length(edges);
num_datasets = length(connectomes);
counts = zeros(numedges-1, num_datasets);

for dataset_idx = 1:num_datasets
    deg(dataset_idx,:) = degrees_und(connectomes{dataset_idx});
    counts(:, dataset_idx) = histcounts(deg(dataset_idx,:), edges);
end

avg_counts = mean(counts,2);
if doPlot
    figure;
    bar(edges(1:end-1), avg_counts, 'FaceColor', [1 .43 .29], 'EdgeColor', [.45 .45 .45]);
    title(sprintf('Average of histograms for %d subjects', num_datasets))
end
end