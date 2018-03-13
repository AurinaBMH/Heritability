% import information
clear all;

cd ('data/general')
filename = 'general_info_restricted.xlsx';
generalInfo = importGeneralInfo(filename);
filename2 = 'general_info_behavioural.xlsx';
[subject,gender] = importGeneralBehav(filename2);

% remove non-twins
nonTwin = find(strcmp(generalInfo.ZygosityGT, ' '));
generalInfo(nonTwin,:) = [];

% find pairs for MZ and DZ twins
% check for family ID; if same, compare mother and father IDs and age (age might differ at testing time).
numSubjects = size(generalInfo,1);
twinTable = array2table(zeros(numSubjects/2,7), 'VariableNames', {'twin1','age1', 'sex1', 'twin2', 'age2', 'sex2', 'zigosity'});
s=1;
for twin1 = 1:numSubjects
    % selected twin should not be already in the table
    % also make sure the first is involved in the calculation
    if  isempty(intersect(generalInfo.Subject(twin1), twinTable.twin1)) ...
            && isempty(intersect(generalInfo.Subject(twin1), twinTable.twin2)) ...
            || twin1 == 1
        is1Twin = generalInfo.ZygosityGT{twin1};
        
        % find same family ID
        familyID = generalInfo.Family_ID{twin1};
        fIND = find(strcmp(familyID, generalInfo.Family_ID));
        % exclude the first twin - there should be only one left
        twin2 = setdiff(fIND,twin1);
        % if only one subject with the same family ID is found compare
        % zigosity info, mother ID, father ID and age
        if length(twin2) == 1
            is2Twin(1) = strcmp(is1Twin, generalInfo.ZygosityGT(twin2));
            is2Twin(2) = generalInfo.Father_ID(twin1)==generalInfo.Father_ID(twin2);
            is2Twin(3) = generalInfo.Mother_ID(twin1)==generalInfo.Mother_ID(twin2);
            is2Twin(4) = generalInfo.Age_in_Yrs(twin1)==generalInfo.Age_in_Yrs(twin2);
            
        else
            % if more than one subject with the same ID is found, print
            % their ID
            fprintf('%d (%d) subject has more than one twin: %s (%s) \n', ...
                generalInfo.Subject(twin1), twin1, generalInfo.Subject(twin2), twin2)
        end
        
        % if all info matches of differs only in age, include them into the
        % table
        if sum(is2Twin) == 4 || (is2Twin(1)==1 && is2Twin(2)==1 && is2Twin(2)==1 && is2Twin(4)==0)
            if is2Twin(1)==1 && is2Twin(2)==1 && is2Twin(2)==1 && is2Twin(4)==0
                % is age differs, print thet on screen
                fprintf('%d (ind %d, age %d) and %d (ind %d, age %d) are twins but age differs\n', ...
                    generalInfo.Subject(twin1), twin1, generalInfo.Age_in_Yrs(twin1), ...
                    generalInfo.Subject(twin2), twin2, generalInfo.Age_in_Yrs(twin2))
            end
            % get subject ID
            twinTable.twin1(s) = generalInfo.Subject(twin1);
            twinTable.twin2(s) = generalInfo.Subject(twin2);
            % find index to get sex info
            ind1 = find(subject == generalInfo.Subject(twin1));
            ind2 = find(subject == generalInfo.Subject(twin2));
            % get age
            twinTable.age1(s) = generalInfo.Age_in_Yrs(twin1);
            twinTable.age2(s) = generalInfo.Age_in_Yrs(twin2);
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
            fprintf('%d (ind %d) subject indo does not match\n', generalInfo.Subject(twin1), twin1)
        end
    end
    
    writetable(twinTable, 'twinCovariates.txt', 'Delimiter', '\t')
    % make variables for MZ twins: ID, age, sex
    MZ_ID(:,1) = twinTable.twin1(twinTable.zigosity==1); 
    MZ_ID(:,2) = twinTable.twin2(twinTable.zigosity==1); 
    
    MZ_age(:,1) = twinTable.age1(twinTable.zigosity==1); 
    MZ_age(:,2) = twinTable.age2(twinTable.zigosity==1); 
    
    MZ_sex(:,1) = twinTable.sex1(twinTable.zigosity==1); 
    MZ_sex(:,2) = twinTable.sex2(twinTable.zigosity==1); 
    
    % make variables for DZ twins: ID, age, sex
    DZ_ID(:,1) = twinTable.twin1(twinTable.zigosity==2); 
    DZ_ID(:,2) = twinTable.twin2(twinTable.zigosity==2); 
    
    DZ_age(:,1) = twinTable.age1(twinTable.zigosity==2); 
    DZ_age(:,2) = twinTable.age2(twinTable.zigosity==2); 
    
    DZ_sex(:,1) = twinTable.sex1(twinTable.zigosity==2); 
    DZ_sex(:,2) = twinTable.sex2(twinTable.zigosity==2); 
    
    % save in matlab file
    save('twinCovariates.mat', 'MZ_ID', 'MZ_age', 'MZ_sex', 'DZ_ID', 'DZ_age', 'DZ_sex'); 

end
cd ../..