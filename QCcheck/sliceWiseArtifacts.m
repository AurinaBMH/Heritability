%% author Aurina 24/07/2017
%% script for detecting slice-wise intensity related artifacts based on Liu et al. (DTIprep)
% read diffusion image
%clear all; close all;
subjects = [107]; % 11 37 43 44 51 105 242];
gradients = [2:11 13:22 24:33 35:44 46:55 57:66]; % select gradients excluding b0.
alpha = 3.5; % threshold for standard deviations
whatCorrelation = 'standard'; %'standard', 'article'
j=1;
n1 = 1;
n2 = 59;
gradRem = zeros(length(subjects),1);
for sub = subjects

    %[d,h] = OpenDiffImage(sub);
    cd(['/gpfs/M2Scratch/Monash076/scratch_shared/GenCog/subjects/1008.2.57.' num2str(sub) '/diffusion/']);
% Read the diffusion images file
%cd(['/gpfs/M2Scratch/Monash076/aurina/Gen_Cog/dwidenoise/1008.2.57.' num2str(subj) ]);
    [h, d] = read('forward.nii'); %dwi_repol_bias_correct forward

    NC = zeros(n2-n1+1,length(gradients));
    NCbad = zeros(n2-n1+1,length(subjects), length(gradients));
    k=1;
    for z=n1:n2
        i=1;
        for grad = gradients

            w1 = squeeze(d(:,:,z,grad));
            w2 = squeeze(d(:,:,z+1,grad));
            switch whatCorrelation
                case 'article'
                    W1 = w1(:);
                    W2 = w2(:);
                    s1 = sum(W1.*W2);
                    s2 = sqrt((sum(W1).*sum(W2)));
                    NC(k,i) = s1/s2;
                case 'standard'
                    NC(k,i) = corr2(w1,w2);
            end
            NCbad(k, j, :) = NC(k,:) < nanmean(NC(k,:))-(alpha.*std(NC(k,:)));
            i=i+1;


        end
        k=k+1;

    end
    matrixBad = squeeze(sum(NCbad,2));
    gradRem(j) = sum(logical(sum(matrixBad,1)));
    figure; imagesc(NC); title (sprintf( 'Subject %d, %d grads to remove', sub, gradRem(j))); xlabel('gradients'); ylabel('slices, z');

    j=j+1;
end
