function [numNodes, leftCortex, leftSubcortex, rightCortex, rightSubcortex] = labelNodes(parcellation)

if strcmp(parcellation, 'aparcANDaseg')
    numNodes = 82;
    leftCortex = 1:34;
    leftSubcortex = 35:41;
    rightCortex = 42:75;
    rightSubcortex = 75:82;
elseif strcmp(parcellation, 'custom100')
    numNodes = 220;
    leftCortex = 1:100;
    leftSubcortex = 101:110;
    rightCortex = 111:210;
    rightSubcortex = 211:220;
elseif strcmp(parcellation, 'custom250')
    numNodes = 530;
    leftCortex = 1:250;
    leftSubcortex = 265;
    rightCortex = 515;
    rightSubcortex = numNodes;
elseif strcmp(parcellation, 'HCPMMP1')
    numNodes = 380;
    leftCortex = 1:180;
    leftSubcortex = 181:190;
    rightCortex = 191:370;
    rightSubcortex = 371:380;
end

end