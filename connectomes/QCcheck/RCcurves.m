function RCcurves(Adj, WhatTypeNetwork,whatNullModel)
if strcmp(WhatTypeNetwork, 'bu')
    Adj = logical(Adj);
end
deg = degrees_und(Adj);
kmax = max(deg);
[PhiNorm,PhiTrue,PhiRand] = RichClubPhiNorm(Adj,kmax,50,100,WhatTypeNetwork,whatNullModel);
%-------------------------------------------------------------------------------
% Compute p values
%-------------------------------------------------------------------------------

pValues = zeros(kmax,1);
for i = 1:kmax
    pValues(i) = mean(PhiTrue(i) <= PhiRand(:,i));
end

isSig = (pValues <= 0.05);
PhiNormMean = zeros(size(PhiTrue));
for i = 1:size(PhiTrue,2)
    PhiNormMean(i) = PhiTrue(i)/nanmean(PhiRand(:,i));
end

figure;
subplot(2,1,2)
h_RCcurve = plot(PhiNormMean, '-','Color','r','LineWidth',2); hold on;
% Significance at p = 0.05 as circles
plot(find(isSig),PhiNormMean(isSig),'o','MarkerEdgeColor','r',...
    'MarkerFaceColor',[1 .41 .38],'LineWidth',1,'MarkerSize',8); xlim([min(deg)-0.5,max(deg)+3]);
xlabel('Degree, k','fontsize',14);
ylabel('\Phi_{norm}', 'fontsize',14);
get(gca, 'XTick');
set(gca, 'FontSize', 12)
box off;
subplot(2,1,1)
histogram(deg, 40, 'FaceColor', [1 .43 .29], 'EdgeColor', [.45 .45 .45]); xlim([min(deg)-0.5,max(deg)+3]);
xlabel('Degree, k','fontsize',14);

box('off')

end