%% Group matrix according to edge length
% combine all length matrices to hhh
% clear all;
% close all;

cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/Final88/');
Parcellation = {'cust100'};
Tract = {'FACT_'};

Name = strcat(Tract{1}, Parcellation{1});
%load (sprintf('%s.mat', Name));

Matrixbin = cell(size(density,2),1);
NumSubj = size(density,2);
BrainPart = {'WithinHemisphere', 'BetweenHemispheres'};

if strcmp(Parcellation, 'aparcANDaseg')

    NumNodes = 82;
    LeftCortex = 34;
    LeftSubcortex = 41;
    RightCortex = 75;
    RightSubcortex = NumNodes;
elseif strcmp(Parcellation, 'cust100')
    NumNodes = 220;
    LeftCortex = 100;
    LeftSubcortex = 110;
    RightCortex = 210;
    RightSubcortex = NumNodes;
elseif strcmp(Parcellation, 'cust250')
    NumNodes = 530;
    LeftCortex = 250;
    LeftSubcortex = 265;
    RightCortex = 515;
    RightSubcortex = NumNodes;

end

matrices = zeros(NumSubj, NumNodes, NumNodes);
% does mean binary density means: the mean of separate densityes or the
% density of mean matrix (when all matrices are averaged before?) or density separately for interhemispheric matrix and intrahemispheric matrix separately.


NumLinks = zeros(NumSubj,1);
NumLinksW = zeros(NumSubj,1);
NumLinksB = zeros(NumSubj,1);

for p = 1:NumSubj

    %Mask = count{p};
    %count{p}(Mask<10) = 0;

    MatrixS1W = count{1,p}(1:LeftSubcortex, 1:LeftSubcortex);
    MatrixS2W = count{1,p}(LeftSubcortex+1:NumNodes, LeftSubcortex+1:NumNodes);
    MatrixW = cat(1, MatrixS1W, MatrixS2W);
    NumLinksW(p) = sum(sum(logical(MatrixW)));

    MatrixS1B = count{1,p}(1:LeftSubcortex, LeftSubcortex+1:NumNodes);
    MatrixS2B = count{1,p}(LeftSubcortex+1:NumNodes, 1:LeftSubcortex);
    MatrixB = cat(1, MatrixS1B, MatrixS2B);
    NumLinksB(p) = sum(sum(logical(MatrixB)));

end
NumOfBinsWithin = 10; %round(sqrt(mean(NumLinksW)));
NumOfBinsBetween = 10; %round(sqrt(mean(NumLinksB)));

MaxNumBins = max([NumOfBinsWithin NumOfBinsBetween]);
KeepLinks = cell(MaxNumBins,4);
% run separately for interhemispheric and intrahemispheric connections.

for run = 1:2
    if run == 1
        NumOfBins = NumOfBinsWithin;
    elseif run == 2
        NumOfBins = NumOfBinsBetween;
    end

    NumOfEdgesS1 = zeros(NumSubj, NumOfBins);
    NumOfEdgesS2 = zeros(NumSubj, NumOfBins);
    MatrixAll1 = zeros(NumSubj, NumOfBins, NumNodes/2, NumNodes/2);
    MatrixAll2 = zeros(NumSubj, NumOfBins, NumNodes/2, NumNodes/2);

    Part = BrainPart(run);

    for i=1:NumSubj
        % for each subject select the appropriate part of the matrix
        switch Part{1}
            case 'WithinHemisphere'
                MatrixS1 = len{1,i}(1:LeftSubcortex, 1:LeftSubcortex);
                %MatrixS1 = tril(MatrixS1);
                MatrixS2 = len{1,i}(LeftSubcortex+1:NumNodes, LeftSubcortex+1:NumNodes);
                %MatrixS2 = tril(MatrixS2);
                Matrix = cat(1, MatrixS1, MatrixS2);

            case 'BetweenHemispheres'
                MatrixS1 = len{1,i}(1:LeftSubcortex, LeftSubcortex+1:NumNodes);
                MatrixS2 = len{1,i}(LeftSubcortex+1:NumNodes, 1:LeftSubcortex);
                Matrix = cat(1, MatrixS1, MatrixS2);

        end
        %define edges for binning
        edges = linspace(min(nonzeros(Matrix(:))),max(Matrix(:)),NumOfBins+1);
        bin = 1:1:size(edges,2);

        % bin length matrix according to edges
        MatrixS1bin = discretize(MatrixS1,edges, 'IncludedEdge','right');
        k = find(isnan(MatrixS1bin))';
        MatrixS1bin(k) = 0;
        MatrixS2bin = discretize(MatrixS2,edges, 'IncludedEdge','right');
        u = find(isnan(MatrixS2bin))';
        MatrixS2bin(u) = 0;


        Matrixbin{i,1} = MatrixS1bin;
        Matrixbin{i,2} = MatrixS2bin;


        % for each bin find how many links fall into each bin for each participant
        for k=1:length(bin)

            [C1] = find(MatrixS1bin == k);
            NumOfEdgesS1(i,k) = size(C1,1);

            [C2] = find(MatrixS2bin == k);
            NumOfEdgesS2(i,k) = size(C2,1);

        end
        % get the average number of edhes in the bin for all subjects
        AvNumOfEdgesS1 = round(mean(NumOfEdgesS1,1));
        AvNumOfEdgesS2 = round(mean(NumOfEdgesS2,1));

    end

    % for each subject and each bin make a "mask" of links within each bin by
    % seting others to 0
    for i=1:NumSubj

        for m=1:length(bin)

            Matrix1 = Matrixbin{i,1};
            Matrix2 = Matrixbin{i,2};

            Matrix1(Matrix1 ~= m) = 0;
            Matrix2(Matrix2 ~= m) = 0;
            MatrixAll1(i,m,:,:) = logical(Matrix1);
            MatrixAll2(i,m,:,:) = logical(Matrix2);

        end
    end



    % find out most common connections in each bin throughout subjects
    for l=1:length(bin)

        switch Part{1}

            case 'BetweenHemispheres'
                if AvNumOfEdgesS1(l) ~=0
                    % for each bin calculate: how popular each connection is
                    PopBin1 = squeeze(mean(MatrixAll1,1));
                    PopMat1 = squeeze(PopBin1(l,:,:));
                    %sort the in decending order to get most common ones
                    [sortedX1,sortingIndices1] = sort((PopMat1(:)),'descend');
                    maxValues1 = sortedX1(1:AvNumOfEdgesS1(l));
                    % set a threshold and randomise indexes for values less than a
                    % threshold
                    Thr1 = min(maxValues1);
                    Test1 = sortedX1(sortedX1>=Thr1);
                    indices1 = randperm(size(Test1,1));
                    %select a required number from randomised indexes
                    indices1 = indices1(1:AvNumOfEdgesS1(l));
                    maxValueIndices1 = sortingIndices1(indices1);

                    % keep the average number of most popular links according to
                    % AvNumOfEdges in a corresponding bin

                    %maxValueIndices1 = sortingIndices1(1:AvNumOfEdgesS1(l));

                    % keep edges with those indexes for each bin
                    [Row1,Col1] = ind2sub(size(Matrix1),maxValueIndices1);
                    BIN1 = zeros(AvNumOfEdgesS1(l),1);
                    BIN1(:) = l;
                    %make a lis of what links in which bin to keep
                    KeepLinks{l,run} = cat(2,BIN1, Row1, Col1);
                    BIN2 = BIN1;
                    run2 = run+1;
                    KeepLinks{l,run2+1} = cat(2,BIN2, Col1, Row1);
                end

            case 'WithinHemisphere'
                if AvNumOfEdgesS1(l) ~=0
                    PopBin1 = squeeze(mean(MatrixAll1,1));
                    PopMat1 = squeeze(PopBin1(l,:,:));
                    %sort the in decending order to get most common ones
                    [sortedX1,sortingIndices1] = sort((PopMat1(:)),'descend');
                    maxValues1 = sortedX1(1:AvNumOfEdgesS1(l));
                    Thr1 = min(maxValues1);
                    Test1 = sortedX1(sortedX1>=Thr1);
                    indices1 = randperm(size(Test1,1));
                    indices1 = indices1(1:AvNumOfEdgesS1(l));
                    maxValueIndices1 = sortingIndices1(indices1);

                    % keep the average number of most popular links according to
                    % AvNumOfEdges in a corresponding bin

                    %maxValueIndices1 = sortingIndices1(1:AvNumOfEdgesS1(l));

                    % keep edges with those indexes for each bin
                    [Row1,Col1] = ind2sub(size(Matrix1),maxValueIndices1);
                    BIN1 = zeros(AvNumOfEdgesS1(l),1);
                    BIN1(:) = l;
                    %make a lis of what links in which bin to keep
                    KeepLinks{l,run} = cat(2,BIN1, Row1, Col1);
                end
                % do the same for the second quarter
                if AvNumOfEdgesS2(l) ~=0
                    PopBin2 = squeeze(mean(MatrixAll2,1));
                    PopMat2 = squeeze(PopBin2(l,:,:));
                    [sortedX2,sortingIndices2] = sort((PopMat2(:)),'descend');
                    %[sortedX2,sortingIndices2] = sort(PopMat2(:),'descend');
                    maxValues2 = sortedX2(1:AvNumOfEdgesS2(l));
                    Thr2 = min(maxValues2);
                    Test2 = sortedX2(sortedX2>=Thr2);
                    indices2 = randperm(size(Test2,1));
                    indices2 = indices2(1:AvNumOfEdgesS2(l));
                    maxValueIndices2 = sortingIndices2(indices2);

                    %maxValueIndices2 = sortingIndices2(1:AvNumOfEdgesS2(l));

                    % keep edges with those indexes for each bin
                    [Row2,Col2] = ind2sub(size(Matrix2),maxValueIndices2);
                    BIN2 = zeros(AvNumOfEdgesS2(l),1);
                    BIN2(:) = l;
                    run2 = run+1;
                    KeepLinks{l,run2+1} = cat(2,BIN2, Row2, Col2);
                end
        end
    end
end



% reconstruct matrices from link files.
% make a mask of what connections to keep combining all bins.
MatrixMasks = cell(4,1);


for side = 1:4
    MatrixSide = zeros(NumNodes/2, NumNodes/2);
    for bine = 1:NumOfBins
        List = KeepLinks{bine, side};
        for r=1:size(List,1)
            % for each bin and link make it 1, if you need to keep a
            % connection. This way a mask for each bin will be made and
            % overlapped

            MatrixSide(List(r,2),List(r,3)) = bine;

        end
    end
    MatrixMasks{side} = MatrixSide;
end

M1 = MatrixMasks{1};
M3 = MatrixMasks{3};
n=size(M1,1);

for b=1:n-1
    for a=b+1:n
        M1(b,a) = M1(a,b);
        M3(b,a) = M3(a,b);
    end
end

MatrixMasks{1} = M1;
MatrixMasks{3} = M3;

Maskup = cat(2,MatrixMasks{1}, MatrixMasks{2});
Maskdw = cat(2, MatrixMasks{4}, MatrixMasks{3});
% combine parts of a mask to make a final mask for the whole matrix
Maskfin = cat(1, Maskup, Maskdw);

Matrix_pop = zeros(NumSubj,NumNodes,NumNodes);

% for each subject keep links according to the mask and make 3D matrix
% with: subject: nodes: nodes
for v=1:NumSubj
    MatrixS = density{1,v};
    Matrix_pop(v,:,:) = MatrixS.*Maskfin;
end

% calculate average weights for those links excluding zeros
Adj = squeeze(sum(Matrix_pop,1) ./ sum(Matrix_pop(:,:,:)~=0));
% replace NaN values with zero
Adj(isnan(Adj)) = 0;

%save([sprintf('%s%s_Lengthgroup', Tract{1}, Parcellation{1})], 'Adj');
