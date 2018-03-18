
function [dens] = plotDensityDistrib(connectomes, doPlot)


% calculate subject-subject correlations
num_datasets = length(connectomes);
dens = zeros(num_datasets,1);

for s1 = 1:num_datasets
    dens(s1) = density_und(connectomes{s1});     
end

if doPlot
    figure; 
    histogram(dens, 100);title('Density distribution')
end
end
