clear all; close all;
%-----------------------------------------------------------------
% Import information about subjects
%-----------------------------------------------------------------
cd ('data/general')
filenameRestricted = 'general_info_restricted.xlsx';
infoALL = importRestricted(filenameRestricted);
filenameGeneral = 'general_info_behavioural.xlsx';
[subject,gender] = importGeneralBehav(filenameGeneral);
% % for all subjects, if they have genetically determined zigosity, use that
% % information, if not, use self-reported information
% for s=1:size(generalInfo,1)
%     if generalInfo.HasGT{s}
%         generalInfo.ZygosityCOMB{s} = generalInfo.ZygosityGT{s};
%     else
%         if strcmp(generalInfo.ZygositySR{s}, 'NotTwin')
%             generalInfo.ZygosityCOMB{s} = ' ';
%         elseif strcmp(generalInfo.ZygositySR{s}, 'NotMZ')
%             generalInfo.ZygosityCOMB{s} = 'DZ';
%         else
%             generalInfo.ZygosityCOMB{s} = generalInfo.ZygositySR{s};
%         end
%     end
% end
%-----------------------------------------------------------------
% Remove subjects that were not genotyped - zygosity is not confirmed
% Create separate variables to store information for twins and non-tins
% Create a twin table to store the info
%-----------------------------------------------------------------
% Find subjects that don't have genetically confirmed zygosity and remove them 
noGT = find(cell2mat(infoALL.HasGT)==0); 
infoALL(noGT,:) = [];

infoTWIN = infoALL; 
infoNONTWIN = infoALL;

% find twins and non twins and creat eseparate viriables for them
NTind = find(strcmp(infoALL.ZygosityGT, ' '));
Tind = find(strcmp(infoALL.ZygosityGT, 'MZ') | strcmp(infoALL.ZygosityGT, 'DZ'));
infoTWIN(NTind,:) = [];
infoNONTWIN(Tind,:) = [];
% find pairs for MZ and DZ twins
% check for family ID; if same, compare mother and father IDs and age (age might differ at testing time).
numSubjects = size(infoTWIN,1);
twinTable = array2table(zeros(numSubjects/2,13), 'VariableNames', {'twin1','age1', 'sex1', 'twin2', 'age2', 'sex2', 'zigosity','sib1','ageS1', 'sexS1','sib2','ageS2', 'sexS2'});
s=1;
for twin1 = 1:numSubjects
    % selected twin should not be already in the table
    % also make sure the first is involved in the calculation
    if  isempty(intersect(infoTWIN.Subject(twin1), twinTable.twin1)) ...
            && isempty(intersect(infoTWIN.Subject(twin1), twinTable.twin2)) ...
            || twin1 == 1
        is1Twin = infoTWIN.ZygosityGT{twin1};
        
        % find same family ID
        familyID = infoTWIN.Family_ID{twin1};
        fIND = find(strcmp(familyID, infoTWIN.Family_ID));
        % exclude the first twin - there should be only one left
        twin2 = setdiff(fIND,twin1);
        % if only one subject with the same family ID is found compare
        % zigosity info, mother ID, father ID and age
        if length(twin2) == 1
            is2Twin(1) = strcmp(is1Twin, infoTWIN.ZygosityGT(twin2));
            is2Twin(2) = infoTWIN.Father_ID(twin1)==infoTWIN.Father_ID(twin2);
            is2Twin(3) = infoTWIN.Mother_ID(twin1)==infoTWIN.Mother_ID(twin2);
            is2Twin(4) = infoTWIN.Age_in_Yrs(twin1)==infoTWIN.Age_in_Yrs(twin2);
            
        else
            % if more than one subject with the same ID is found, print
            % their ID
            fprintf('%d (%d) subject has more than one twin: %s (%s) \n', ...
                infoTWIN.Subject(twin1), twin1, infoTWIN.Subject(twin2), twin2)
        end
        
        % if all info matches of differs only in age, include them into the
        % table
        if sum(is2Twin) == 4 || (is2Twin(1)==1 && is2Twin(2)==1 && is2Twin(2)==1 && is2Twin(4)==0)
            if is2Twin(1)==1 && is2Twin(2)==1 && is2Twin(2)==1 && is2Twin(4)==0
                % is age differs, print thet on screen
                fprintf('%d (ind %d, age %d) and %d (ind %d, age %d) are twins but age differs\n', ...
                    infoTWIN.Subject(twin1), twin1, infoTWIN.Age_in_Yrs(twin1), ...
                    infoTWIN.Subject(twin2), twin2, infoTWIN.Age_in_Yrs(twin2))
            end
            % get subject ID
            twinTable.twin1(s) = infoTWIN.Subject(twin1);
            twinTable.twin2(s) = infoTWIN.Subject(twin2);
            % find index to get sex info
            ind1 = find(subject == infoTWIN.Subject(twin1));
            ind2 = find(subject == infoTWIN.Subject(twin2));
            % get age
            twinTable.age1(s) = infoTWIN.Age_in_Yrs(twin1);
            twinTable.age2(s) = infoTWIN.Age_in_Yrs(twin2);
            % get sex
            if strcmp(gender(ind1), 'M')
                twinTable.sex1(s) = 1;
            elseif strcmp(gender(ind1), 'F')
                twinTable.sex1(s) = 2;
            end
            
            if strcmp(gender(ind2), 'M')
                twinTable.sex2(s) = 1;
            elseif strcmp(gender(ind2), 'F')
                twinTable.sex2(s) = 2;
            end
            
            % give 1 label for MZ
            if strcmp(is1Twin, 'MZ')
                twinTable.zigosity(s) = 1;
            % give 2 label for DZ
            elseif strcmp(is1Twin, 'DZ')
                twinTable.zigosity(s) = 2;
            end
            s=s+1;
            
        else
            % if more than age missmatches, print on screen
            fprintf('%d (ind %d) subject indo does not match\n', infoTWIN.Subject(twin1), twin1)
        end
    end

end
%writetable(twinTable, 'twinCovariates.txt', 'Delimiter', '\t')
%-----------------------------------------------------------------
% For each pair of confirmed twins, find if there are non-twin siblings and
% add them into the table
%-----------------------------------------------------------------

sib=1;
cell_array = infoNONTWIN.Family_ID;
for t=1:size(twinTable,1)
    % for each twin from the pairs, find family ID in the intial data
    ind1 = find(infoTWIN.Subject==twinTable.twin1(t));
    ind2 = find(infoTWIN.Subject==twinTable.twin2(t));
    % compare if they match - they should by default - we confirmed that
    % already, but just in case
    if strcmp(infoTWIN.Family_ID{ind1},infoTWIN.Family_ID{ind2})
        cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
        string = infoTWIN.Family_ID{ind1};
        % find if there are other subjects that have the same family ID
        siblingIND = find(cellfun(cellfind(string),cell_array));
        % if there is one additional signling add it as S1 and give S2 nan
        % values
        if length(siblingIND)==1
            
            twinTable.sib1(t) = infoNONTWIN.Subject(siblingIND);
            twinTable.ageS1(t) = infoNONTWIN.Age_in_Yrs(siblingIND);
            twinTable.sib2(t) = NaN;
            twinTable.ageS2(t) = NaN;
            twinTable.sexS2(t) = NaN; 
            % find index to get sex info
            indS = find(subject == infoNONTWIN.Subject(siblingIND));
            if strcmp(gender(indS), 'M')
                twinTable.sexS1(t) = 1;
            elseif strcmp(gender(indS), 'F')
                twinTable.sexS1(t) = 2;
            end
            sib=sib+1;
        % if there are more than 1, give values to both
        elseif length(siblingIND)>1
            fprintf('Twins in row %d have %d siblings\n', t, length(siblingIND))
            % choose one sibling
            twinTable.sib1(t) = infoNONTWIN.Subject(siblingIND(1));
            twinTable.ageS1(t) = infoNONTWIN.Age_in_Yrs(siblingIND(1));
            twinTable.sib2(t) = infoNONTWIN.Subject(siblingIND(2));
            twinTable.ageS2(t) = infoNONTWIN.Age_in_Yrs(siblingIND(2));
            % find index to get sex info
            indS1 = find(subject == infoNONTWIN.Subject(siblingIND(1)));
            indS2 = find(subject == infoNONTWIN.Subject(siblingIND(2)));
            if strcmp(gender(indS1), 'M')
                twinTable.sexS1(t) = 1;
            elseif strcmp(gender(indS1), 'F')
                twinTable.sexS1(t) = 2;
            end
            
            if strcmp(gender(indS2), 'M')
                twinTable.sexS2(t) = 1;
            elseif strcmp(gender(indS2), 'F')
                twinTable.sexS2(t) = 2;
            end
            
            sib=sib+1;
        else
            twinTable.sib1(t) = NaN;
            twinTable.ageS1(t) = NaN;
            twinTable.sexS1(t) = NaN; 
            twinTable.sib2(t) = NaN;
            twinTable.ageS2(t) = NaN;
            twinTable.sexS2(t) = NaN; 
        end
    else
        warning(sprintf('Family IDs do not match for twins in row %d\n', t))
    end
end
fprintf('Total number of twins with a sigling is %d\n', sib)


% make variables for MZ twins and their siblings: ID, age, sex
MZ_ID(:,1) = twinTable.twin1(twinTable.zigosity==1);
MZ_ID(:,2) = twinTable.twin2(twinTable.zigosity==1);
MZ_ID(:,3) = twinTable.sib1(twinTable.zigosity==1);
MZ_ID(:,4) = twinTable.sib2(twinTable.zigosity==1);

MZ_age(:,1) = twinTable.age1(twinTable.zigosity==1);
MZ_age(:,2) = twinTable.age2(twinTable.zigosity==1);
MZ_age(:,3) = twinTable.ageS1(twinTable.zigosity==1);
MZ_age(:,4) = twinTable.ageS2(twinTable.zigosity==1);

MZ_sex(:,1) = twinTable.sex1(twinTable.zigosity==1);
MZ_sex(:,2) = twinTable.sex2(twinTable.zigosity==1);
MZ_sex(:,3) = twinTable.sexS1(twinTable.zigosity==1);
MZ_sex(:,4) = twinTable.sexS2(twinTable.zigosity==1);

% make variables for DZ twins and their siblings: ID, age, sex
DZ_ID(:,1) = twinTable.twin1(twinTable.zigosity==2);
DZ_ID(:,2) = twinTable.twin2(twinTable.zigosity==2);
DZ_ID(:,3) = twinTable.sib1(twinTable.zigosity==2);
DZ_ID(:,4) = twinTable.sib2(twinTable.zigosity==2);


DZ_age(:,1) = twinTable.age1(twinTable.zigosity==2);
DZ_age(:,2) = twinTable.age2(twinTable.zigosity==2);
DZ_age(:,3) = twinTable.ageS1(twinTable.zigosity==2);
DZ_age(:,4) = twinTable.ageS2(twinTable.zigosity==2);


DZ_sex(:,1) = twinTable.sex1(twinTable.zigosity==2);
DZ_sex(:,2) = twinTable.sex2(twinTable.zigosity==2);
DZ_sex(:,3) = twinTable.sexS1(twinTable.zigosity==2);
DZ_sex(:,4) = twinTable.sexS2(twinTable.zigosity==2);

% save in matlab file
save('twinCovariates.mat', 'MZ_ID', 'MZ_age', 'MZ_sex', 'DZ_ID', 'DZ_age', 'DZ_sex');

%-----------------------------------------------------------------
% To remove twins that don't have DWI, load any connectome file
%-----------------------------------------------------------------
whatDWI = 'HCP';
parc = 'HCP';
tract = 'iFOD2';
weight1 = 'standard';
brainPart = 'wholeBrain'; 
strRem = 10;

[A, matrices, coordinates, avWeight, SUBjects] = giveConnDATA(whatDWI,parc,tract,weight1,brainPart,strRem); 
% remove subjects that don't have DWI. If one in a pair doesn't have it,
% remove both
i=1;
for MZ=1:size(MZ_ID,1)
    t1 = intersect(SUBjects, MZ_ID(MZ,1));
    t2 = intersect(SUBjects, MZ_ID(MZ,2));
    t3 = intersect(SUBjects, MZ_ID(MZ,3));
    t4 = intersect(SUBjects, MZ_ID(MZ,4));
    
    if isempty(t1) || isempty(t2)
        remMZ(i) = MZ;
        i=i+1;
    end
    % if there is no additional sibling,replace their values with NaNs
    if isempty(t3)
        MZ_ID(MZ,3) = NaN; 
        MZ_age(MZ,3) = NaN; 
        MZ_sex(MZ,3) = NaN; 
    end
    
    if isempty(t4)
        MZ_ID(MZ,4) = NaN; 
        MZ_age(MZ,4) = NaN; 
        MZ_sex(MZ,4) = NaN; 
    end
end
% do same for DZ
j=1;
for DZ=1:size(DZ_ID,1)
    t1 = intersect(SUBjects, DZ_ID(DZ,1));
    t2 = intersect(SUBjects, DZ_ID(DZ,2));
    t3 = intersect(SUBjects, DZ_ID(DZ,3));
    t4 = intersect(SUBjects, DZ_ID(DZ,4));
    
    if isempty(t1) || isempty(t2)
        remDZ(j) = DZ;
        j=j+1;
    end 
    % if there is no additional sibling,replace their values with NaNs
    if isempty(t3)
        DZ_ID(DZ,3) = NaN; 
        DZ_age(DZ,3) = NaN; 
        DZ_sex(DZ,3) = NaN; 
    end
    if isempty(t4)
        DZ_ID(DZ,4) = NaN;
        DZ_age(DZ,4) = NaN; 
        DZ_sex(DZ,4) = NaN; 
    end
end

% remove pairs if one of the twins is missing diffusion data
MZ_ID(remMZ,:) = [];
DZ_ID(remDZ,:) = [];

MZ_age(remMZ,:) = [];
DZ_age(remDZ,:) = [];

MZ_sex(remMZ,:) = [];
DZ_sex(remDZ,:) = [];

fprintf('MZ: %d out of %d pairs have 1 additional sibling\n', length(find(~isnan(MZ_ID(:,3)))),size(MZ_ID,1));
fprintf('MZ: %d out of %d pairs have 2 additional siblings\n', length(find(~isnan(MZ_ID(:,4)))),size(MZ_ID,1));

fprintf('DZ: %d out of %d pairs have 1 additional sibling\n', length(find(~isnan(DZ_ID(:,3)))),size(DZ_ID,1));
fprintf('DZ: %d out of %d pairs have 2 additional siblings\n', length(find(~isnan(DZ_ID(:,4)))),size(DZ_ID,1));

save('twinCovariatesDWI.mat', 'MZ_ID', 'MZ_age', 'MZ_sex', 'DZ_ID', 'DZ_age', 'DZ_sex');
cd ../..

