% this script reorders random parcellation according to nodes in anatomical
% parcellation

function Matrix = fRandom2AnatomicalParcellations_annot(Parcellation, Subcparc, Adj)

cd ('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Connectomes/Labels/');
% choose the parcellation
%Parcellation = {'custom200'};
Sides = {'lh', 'rh'};
hem=1;
order = cell(length(Sides),1);

% choose tractography methos (iFOD2 or FACT)
%Method = {'FACT'};
% choose subcortical parcellation
%Subcparc = {'ANDaseg20'};
% define number of nodes in subcortex
if strcmp(Subcparc, 'ANDaseg20')
    NumSubcNodes = 10;
elseif strcmp(Subcparc, 'ANDaseg')
    NumSubcNodes = 7;
elseif strcmp(Subcparc, 'ANDfslatlas')
    NumSubcNodes = 14;
elseif strcmp(Subcparc, 'ANDfslatlas20')
    NumSubcNodes = 10;
elseif strcmp(Subcparc, 'ANDaseg30')
    NumSubcNodes = 15;
elseif strcmp(Subcparc, 'NO')
    NumSubcNodes = 0;
end

% choose weight (count or density)
%Weight = {'count'};


for side = Sides
    % load anotation files separately for right and left cortex
    AnatFileName = sprintf('%s.aparc.annot', side{1});
    RandFileName = sprintf('%s.%s.annot', side{1}, Parcellation);
    
    [aparcvertices, aparclabel, aparccolortable] = read_annotation(AnatFileName);
    [randvertices, randlabel, randcolortable] = read_annotation(RandFileName);
    
    
    anat = aparccolortable.table([2:4 6:end],5); % exclude unknown and corpuscallosum, results in 34 nodes (for 1 hemisphere)
    rand = randcolortable.table(2:end,5); % exclude unknown
    
    NumNodesAnat = length(anat);
    NumNodesRand = length(rand);
    NumNodes = NumNodesRand*2 + 2*NumSubcNodes;
    
    
    r1 = cell(NumNodesRand,1);
    for i=1:size(rand,1)
        noderand = rand(i);
        r1{i} = find(randlabel == noderand);
    end
    
    
    r2 = cell(NumNodesAnat,1);
    for j=1:size(anat,1)
        nodeanat = anat(j);
        r2{j} = find(aparclabel == nodeanat);
    end
    
    % find the intersection of anatomical and random lists
    inter = zeros(NumNodesRand,NumNodesAnat);
    for k=1:NumNodesRand
        randList = r1{k};
        for l=1:NumNodesAnat
            Anatlist = r2{l};
            inter(k,l) = length(intersect(randList,Anatlist));
        end
    end
    
    % for each random node find a corresponding anatomical node according
    % to the max intersection value.
    [M,I] = max(inter,[],2);
    randList = 1:1:NumNodesRand;
    list = [randList', I];
    order{hem} = sortrows(list,2);
    hem=hem+1;
end
% go to connectomes
% cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Connectomes/10_50_mlnStreamlines/');
% % load relevant connectome as defined in the beginning of the script
% load(sprintf('%s_%s%s_CommonConnections_%s.mat', Method{1}, Parcellation{1}, Subcparc{1}, Weight{1}));
% %choose group matrix (4 - connections present in 60% subjects)
%Adj = squeeze(Adj_all(4,:,:));

%cortex-cortex (take a part of a matrix and reorder)
orderL = order{1};
orderR = order{2};

%select top rows of the matrix and reorder according to left
Lhoz1 = Adj(1:NumNodesRand, 1:end);
Lhoz1 = Lhoz1(orderL(:,1),1:end);
Lhoz2 = Adj(NumNodesRand+1:NumNodesRand+NumSubcNodes, 1:end);

L = vertcat(Lhoz1, Lhoz2);

% select bottom rows of a matrix and reorder according to right
Rhoz1 = Adj(NumNodesRand+NumSubcNodes+1:end, 1:end);
Rhoz2 = Adj(NumNodes-NumSubcNodes+1:end, 1:end);
Rhoz1 = Rhoz1(orderR(:,1), 1:end);
R = vertcat(Rhoz1, Rhoz2);

% combine to a matrix
RL = vertcat(L,R);

%select vertical strips which do not need reordering
vert1 = RL(1:end, NumNodesRand+1: NumNodesRand+NumSubcNodes);
vert2 = RL(1:end, NumNodes-NumSubcNodes+1:end);

% select and reorder 1st column according to left
Lvert = RL(1:end, 1:NumNodesRand);
Lvert = Lvert(1:end, orderL(:,1));

% select and reorder 2nd column according to right
Rvert = RL(1:end, NumNodesRand+NumSubcNodes+1:NumNodes-NumSubcNodes);
Rvert = Rvert(1:end, orderR(:,1));

% combine reordered parts with vert1 and vert2 strips
Matrix = horzcat(Lvert, vert1, Rvert, vert2);
end




