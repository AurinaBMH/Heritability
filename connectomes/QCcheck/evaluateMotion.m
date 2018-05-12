% arrange data
clear all; close all; 

motion = CopyHCP_motion(); 

scatterColor = [.98 .85 .37]; 
histoColor = [.28 .75 .57]; 
edgeColor = [0.45 0.45 0.45]; 

% plot SNR scatterplots for all combinations
figure; set(gcf,'color','w');
subplot(2,3,1); scatter(motion.SNR(:,1), motion.SNR(:,2), 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor);xlabel('b0'); ylabel('b1000');
subplot(2,3,2); scatter(motion.SNR(:,1), motion.SNR(:,3), 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor);xlabel('b0'); ylabel('b2000');
subplot(2,3,3); scatter(motion.SNR(:,1), motion.SNR(:,4), 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor);xlabel('b0'); ylabel('b3000');
subplot(2,3,4); scatter(motion.SNR(:,2), motion.SNR(:,3), 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor);xlabel('b1000'); ylabel('b2000');
subplot(2,3,5); scatter(motion.SNR(:,2), motion.SNR(:,4), 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor);xlabel('b1000'); ylabel('b3000');
subplot(2,3,6); scatter(motion.SNR(:,3), motion.SNR(:,4), 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor);xlabel('b2000'); ylabel('b3000');

% make variables
m_FD = zeros(length(motion.SUBS),1); 
m_rot = zeros(length(motion.SUBS),1); 
m_tran = zeros(length(motion.SUBS),1); 
FDcheck = ones(length(motion.SUBS),1); 

% check if there are subjects that have any value of FD exceeding 2
for i=1:length(motion.SUBS)
    
    % check if there are subjects that have any value of FD exceeding 2
    if sum(motion.LR_FD(i,:)>2)>1 || sum(motion.RL_FD(i,:)>2)>1
        FDcheck(i) = 0;
    end
    
    % calculate mean FD, rotation and translation across left and right
    m_FD(i) = mean(motion.LR_FD(i,:)) + mean(motion.RL_FD(i,:))/2; 
    
   
    m_rot(i) = mean(motion.LR_rotation(i,:)) + mean(motion.RL_rotation(i,:))/2; 
    m_tran(i) = mean(motion.LR_translation(i,:)) + mean(motion.RL_translation(i,:))/2; 
    % calculate the mean of SNR for all b values (they are on different
    % scales but correlated) - doesn't really make sense to average this
    % way
   % m_SNR(i) = mean(motion.SNR(i,:),2); 
    

end

for i=1:18
m(i) = mean(motion.LR_FD(:,i)) + mean(motion.RL_FD(:,i))/2;
s(i) = std(motion.LR_FD(:,i)) + std(motion.RL_FD(:,i))/2;
end

figure; set(gcf,'color','w');errorbar(m,s, 'Color', [.45 .4 .45]); 
hold on; scatter(1:18,m, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', histoColor); 
xlabel('time-points (b0)'); ylabel('mean FD')
box off; 

% label subjects that have a measure exceeding 2 SD
ROTcheck = ~(m_rot> mean(m_rot) + 2*std(m_rot)); 
TRAcheck = ~(m_tran> mean(m_tran) + 2*std(m_tran));
FD2check = ~(m_FD> mean(m_FD) + 2*std(m_FD));
%SNRcheck = ~(m_SNR> mean(m_SNR) + 2*std(m_SNR));

% Evaluate SNR for each b value separately
SNRcheckSep = ones(length(motion.SUBS),4); 
for b=1:4
SNRcheckSep(:,b) = ~(motion.SNR(:,b)< mean(motion.SNR(:,b)) - 2*std(motion.SNR(:,b)));
end

% plot all those measures 
[proportionOutliers,eddyOutput, outliers, eddyMotion] = CopyHCPeddy(); 
allmotion = horzcat(SNRcheckSep, ROTcheck, TRAcheck, FD2check, FDcheck, ~eddyMotion); 
nrSubjects = length(motion.SUBS) - sum(allmotion,1); 
figure; imagesc(allmotion); colormap(hot); caxis([-1 2])
set(gcf,'color','w');
xticks([1:9])
xticklabels({sprintf('SNRb0 (%d)', nrSubjects(1)), sprintf('SNRb1000 (%d)', nrSubjects(2)), ...
    sprintf('SNRb2000 (%d)', nrSubjects(3)), sprintf('SNRb3000 (%d)', nrSubjects(4)), ...
    sprintf( 'rotation (%d)', nrSubjects(5)), sprintf('translation (%d)', nrSubjects(6)), ...
    sprintf('meanFD (%d)',nrSubjects(7)), sprintf('any FD > 2 (%d)', nrSubjects(8)), ...
    sprintf('eddy outliers(%d)', nrSubjects(9))})
xtickangle(45)
ylabel('Subjects')

% plot them
%figure; 
%subplot(2,2,1); histogram(FDcheck, 20); title('Any FD value > 2');
%subplot(2,2,2); histogram(FD2check, 20); title('mean FD value > 2SD');
%subplot(2,2,3); histogram(ROTcheck, 20); title('mean rotation value > 2SD');
%subplot(2,2,4); histogram(TRAcheck, 20); title('mean translation value > 2SD');

% are same subjects high on all motion parameters
motionALL = FD2check+ROTcheck+TRAcheck; 
figure; histogram(motionALL, 'EdgeColor', edgeColor, 'FaceColor', histoColor); 
set(gcf,'color','w');
title('Number of subjects passing QC criteria on motion (FD, translation, rotation)')
xticks([0 1 2 3])
xticklabels({sprintf('Fail all 3 (%d)', length(find(motionALL==0))), ... 
    sprintf('Fail on 2 (%d)',length(find(motionALL==1))), ... 
    sprintf('Fail on 1 (%d)',length(find(motionALL==2))), ...
sprintf('Pass on all (%d)', length(find(motionALL==3)))})

% plot histograms for motion measures
figure; set(gcf,'color','w');
subplot(1,3,1); histogram(m_FD, 50, 'EdgeColor', edgeColor, 'FaceColor', histoColor); title ('mean FD')
subplot(1,3,2); histogram(m_rot, 50, 'EdgeColor', edgeColor, 'FaceColor', histoColor); title ('mean rotation')
subplot(1,3,3); histogram(m_tran, 50,'EdgeColor', edgeColor, 'FaceColor', histoColor); title ('mean translation')

figure; set(gcf,'color','w');
for i=1:4
subplot(2,2,i); histogram(motion.SNR(:,i), 'EdgeColor', edgeColor, 'FaceColor', histoColor); 
title(sprintf('SNR at b %d', motion.b_vals(i)))
end

% check if rotation and translation are correlated
figure; 
subplot(1,3,1); scatter(m_rot, m_tran, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('Rotation'); ylabel('Translation'); 
[r_rt,p_rt] = corr(m_rot, m_tran, 'type', 'Spearman'); 
title(sprintf('r = %s', num2str(r_rt)))

subplot(1,3,2); scatter(m_FD, m_tran, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('FD'); ylabel('Translation'); 
[r_FDt,p_FDt] = corr(m_FD, m_tran, 'type', 'Spearman'); 
title(sprintf('r = %s', num2str(r_FDt)))

subplot(1,3,3); scatter(m_FD, m_rot, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('FD'); ylabel('Rotation'); 
[r_FDr,p_FDr] = corr(m_FD, m_rot, 'type', 'Spearman'); 
title(sprintf('r = %s', num2str(r_FDr)))

% correlation between hean motion and SNR

figure;set(gcf,'color','w');
subplot(3,2,1); scatter(motion.SNR(:,1), m_FD, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('SNR b0'); ylabel('mean FD'); 
[rval, pval] = corr(motion.SNR(:,1), m_FD, 'type', 'Spearman'); 
str = sprintf('r = %s, p=%s',num2str(rval), num2str(pval)); 
title(str); 

subplot(3,2,2); scatter(motion.SNR(:,4), m_FD, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('SNR b3000'); ylabel('mean FD'); 
[rval, pval] = corr(motion.SNR(:,4), m_FD, 'type', 'Spearman'); 
str = sprintf('r = %s, p=%s',num2str(rval), num2str(pval)); 
title(str);

subplot(3,2,3); scatter(motion.SNR(:,1), m_tran, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('SNR b0'); ylabel('mean translation'); 
[rval, pval] = corr(motion.SNR(:,1), m_tran, 'type', 'Spearman'); 
str = sprintf('r = %s, p=%s',num2str(rval), num2str(pval)); 
title(str);

subplot(3,2,4); scatter(motion.SNR(:,4), m_tran, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('SNR b3000'); ylabel('mean translation'); 
[rval, pval] = corr(motion.SNR(:,4), m_tran, 'type', 'Spearman'); 
str = sprintf('r = %s, p=%s',num2str(rval), num2str(pval));
title(str);

subplot(3,2,5); scatter(motion.SNR(:,1), m_rot, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('SNR b0'); ylabel('mean rotation'); 
[rval, pval] = corr(motion.SNR(:,1), m_rot, 'type', 'Spearman'); 
str = sprintf('r = %s, p=%s',num2str(rval), num2str(pval));
title(str);

subplot(3,2,6); scatter(motion.SNR(:,4), m_rot, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('SNR b3000'); ylabel('mean rotation'); 
[rval, pval] = corr(motion.SNR(:,4), m_rot, 'type', 'Spearman'); 
str = sprintf('r = %s, p=%s',num2str(rval), num2str(pval)); 
title(str);

% compare eddy motion to FD
M = load('motion_parameters.mat');
% calculate the average on b0 volumes
load('Volume_index.mat'); 
everySecond = true;  
b0 = [1 17 33 49 65 81 96 113 129 145 161 177 192 209 225 241 257 273];
[b0IND, a] = find(volume_ind==b0);
motionEddy = M.motion;

numSubj = length(motion.SUBS);

% calculate the average of motion
avMotion = zeros(numSubj,2);
avMotionb0 = zeros(numSubj,1);
secondIND = 1:2:576;
secondINDb0 = 1:2:36;
secondINDb01 = 2:2:36;

for s=1:numSubj
    
    b0motion = motionEddy{s}(b0IND,2);
    
    if everySecond
        allmotion = motionEddy{s}(secondIND,2);
        b0motion = (mean(b0motion(secondINDb0)) + mean(b0motion(secondINDb01)))/2;
    else
        allmotion = motionEddy{s}(:,2);
    end
    avMotionb0(s) = b0motion;
    avMotion(s,:) = mean(allmotion,1);
   
end

b0motion = avMotion(b0IND,1);

figure; 
subplot(1,3,1); scatter(m_FD, avMotionb0, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('FD'); ylabel('Eddy motion')
subplot(1,3,2); scatter(m_tran, avMotionb0, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('Translation'); ylabel('Eddy motion')
subplot(1,3,3); scatter(m_rot, avMotionb0, 'MarkerEdgeColor', edgeColor, 'MarkerFaceColor', scatterColor); xlabel('Rotation'); ylabel('Eddy motion')

% plot the proportion of 





