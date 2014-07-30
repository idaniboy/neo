
% your matlab & clips directory
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
nClips=length(clipDirNames);

% test subject
ii = 12;

% init vars for test subject
clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
nIClips=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));
nIIClips=length(dir(strcat(clipLoc,filesep,'*_interictal_segment_*.mat')));
nIIIClips = nIClips + nIIClips;
lfS = clipScores{ii}(1,:);
sS = clipScores{ii}(2,:);

% Create true labels vector
labels = [ones(nIClips,1);zeros((nIIClips),1)];    

% Compute ROC for SVM
posclass = 1;% Define positive class
[X,Y,T,AUC] = perfcurve(labels,lfS,posclass);% Generate AUC

% Plot ROC Curve
figure
plot(X,Y)
xlabel('False positive rate'); ylabel('True positive rate')
title(strcat('ROC, Subj_',num2str(ii)))
annotation('textbox',...
    [0.9 0.1 0.1 0.1],...
    'String',{'AUC:',num2str(AUC)},...
    'FontSize',12,...
    'FontName','Arial',...
    'LineStyle','--',...
    'EdgeColor',[0 0 0],...
    'LineWidth',1,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);
axis square; set(gcf,'color','w');drawnow;    
save(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(ii)),'cSvmM');
drawnow;


% Generate threshold based off of AUC
[~,~,~,~,OPTROCPT] = perfcurve(labels,lfS,posclass);
optTPR = OPTROCPT(1,2); %optimal TPR
[TPR,~,T,~] = perfcurve(labels,lfS,posclass,'Xcrit','TPR');
ind = find(TPR==optTPR);
searchScore = mean(T(ind));
disp(strcat('searchScore for Patient8:',num2str(searchScore)));