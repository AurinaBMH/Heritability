% make a txt file for an edge
clear all; close all;
parcellation = 'HCPMMP1'; % 'HCPMMP1' , 'custom200';
tract = 'iFOD2';
sift = 'SIFT2';
weight = 'FA'; % 'FA', 'standard'

cd ('data/general')

load('twinCovariatesDWI.mat')
cd ..
cd ('connectomes')

Conn = load(sprintf('%sANDfslatlas20_acpc_%s_%s_%s_structnets.mat', parcellation, tract, sift, weight));
SUBS = Conn.SUBS; 
connectomes = Conn.ADJS;
numSubj = length(SUBS); 
edgeMatr = []; 
for s=1:numSubj 
    m = maskuHalf(connectomes{s});
    edgeMatr(:,s) = m(~isnan(m)); 
end

numEdges = size(edgeMatr,1); 
Output_MZ = zeros(size(MZ_ID,1),2, numEdges); 

for i=1:size(MZ_ID,1)
        ind1 = find(SUBS==MZ_ID(i,1)); 
        ind2 = find(SUBS==MZ_ID(i,2)); 
    for edg = 1:numEdges
    Output_MZ(i,1,edg) = edgeMatr(edg,ind1); 
    Output_MZ(i,2,edg) = edgeMatr(edg,ind2); 
    end
end

Output_DZ = zeros(size(DZ_ID,1),2, numEdges); 
for i=1:size(DZ_ID,1)
        ind1 = find(SUBS==DZ_ID(i,1)); 
        ind2 = find(SUBS==DZ_ID(i,2)); 
    for edg = 1:numEdges
    Output_DZ(i,1,edg) = edgeMatr(edg,ind1); 
    Output_DZ(i,2,edg) = edgeMatr(edg,ind2); 
    end
end

cd ..
cd ('output')
fileName = sprintf('twinEdges_%s_%s_%s.mat', parcellation, tract, weight); 
save(fileName, 'Output_MZ', 'Output_DZ');  