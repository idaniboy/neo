%% creates feature vectors from low frequency features of seizure data
% derived from lfSDforClips
% chooses feature based off performance on training data
% plots AUC curves for each subject
% non-normalized LFS on row 1
% normalizes samples on 0 to 1 scale for row 2

%your matlab directory
%matdir='/Users/idaniboy/Documents/MATLAB/';
%matdir='C:\Users\paul\Documents\MATLAB\';
%clips dir(relative to matlab, leave off^ that slash, as below. these get mushed together)
%clipsdir='kaggleShazam\clips\';
%clipsdir='kaggleClips/ictalSegs/';
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';

%requires data to be in separate folders per patient with _ in folder name
a=dir(strcat(matdir,clipsdir,'*_*'));
%clipDirNames={'Dog_1'};%{a.name};
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;

% clipDirNames={'Dog_2'};%,'Dog_2','Dog_3'};
nClips=length(clipDirNames);

% feature vector order sequence; in future automate discovery of this
% using lfSDforClipsPreTesting
fvOrder = {'mll','smmSco','rmsm','tmmSco','mll','mll','mspSco','rmsm','mll','rmsm','rmsm','rmsm'};

% initialize clipScoresrray to score lfScores
clipScores = cell(12,1);
    
%% find lfsScore,
tic
s = warning('off','all'); % turn all warnings off
if true,
    disp('Warnings off while running this code!') %lots of warnings from peakPoints not finding a peak
    warning('Here is a warning that doesn''t have an id.')
end

for nthDir=7:7;  
    %load concatenated_ filenames
    clipLoc=strcat(matdir,clipsdir,clipDirNames{nthDir});
    nIClips=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));
    nIIClips=length(dir(strcat(clipLoc,filesep,'*_interictal_segment_*.mat')));
    nTSClips=length(dir(strcat(clipLoc,filesep,'*_test_segment_*.mat')));

    %initialize parameters
    [~,~,freq,channels]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_segment_',num2str(1),'.mat'));
    freq = ceil(freq);
    numCh = length(fieldnames(channels));
    data = zeros(numCh,ceil(freq));

    %choose method (one per subject)
    %scoMeth = fvOrder{nthDir}; 

    %if multiple methods per subject
    features = {'mspSco'};

    %data store for lfS
    %lfS = zeros(nIClips+nIIClips,length(features));
    %data stores
    TSlfS = zeros(nTSClips,length(features)); %interictal low freq score
        
    for nthMeth=1:length(features);


        disp(['Finding ' features{nthMeth} ' searchScore for: ' clipLoc])

        % Scores for interictal clips
        for i = 1:nTSClips;
            [data,~,freq,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_test_segment_',num2str(i),'.mat'));
            TSlfS(i,nthMeth) = lfScoreGet(data,features{nthMeth},numCh,freq,clipLoc);
        end
        
        % Sort scores ascending order
        for j=1:length(features)
            sortTSlfS(:,j)=sort(TSlfS(:,j));
            figure;plot(sortTSlfS(:,j));title(strcat(clipDirNames{nthDir},features(nthMeth)));axis square; set(gcf,'color','w');drawnow;
        end
    end
end

break

%% SVM
sS = (1-testSS(7,1:size(lfS,1),1))';
% Scores vector == lfS
%X = [lfS(:,1),lfS(:,2),sS];
X = [lfS(:,1),sS];

% Create true labels vector
labels = [ones(nIClips,1);zeros((nIIClips),1)];    

% train svm model
SVMModel = fitcsvm(X,labels,'ClassNames',[0 1],'KernelFunction','rbf','KernelScale','auto');
%cSvmM = crossval(SVMModel);
cSvmM = compact(SVMModel);
cSvmM = fitSVMPosterior(cSvmM,X,labels); %returns probabilities

% Save Model (overwrites previous model!)
save(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(nthDir)),'cSvmM');

% predict w/SVM
[predLabel,postProbs] = predict(cSvmM,X);
qS = postProbs(:,2);    

% Plot results
figure('Position', [100, 500, 1100, 400]);
subplot(1,3,1);plot(qS);title(strcat('SVM Post-Prob. Subj_',clipDirNames(nthDir)));axis square;axis tight;set(gcf,'color','w');drawnow
subplot(1,3,2);
for i=1:size(predLabel,1)
    if predLabel(i) ~= labels(i)
        colr = 'r';
    else
        if labels(i) == 1
            colr = 'g';
        else
            colr = 'b';
        end
    end
    scatter3(sS(i),lfS(i,1),lfS(i,2),colr);xlabel('searchScore');ylabel('slopeScore'),zlabel('phaseScore');title(strcat('SVM Results. Subj_',clipDirNames(nthDir)));hold on
end
axis square;set(gcf,'color','w');  drawnow;  


% Plot training ROC curve
% Create true labels vector
labels = [ones(nIClips,1);zeros((nIIClips),1)];    

% Compute ROC for SVM
posclass = 1;% Define positive class
[X,Y,T,AUC] = perfcurve(labels,qS,posclass);% Generate AUC

% Plot ROC Curve
subplot(1,3,3)

plot(X,Y)
xlabel('False positive rate'); ylabel('True positive rate')
title(strcat('ROC, Subj_',num2str(nthDir)))
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
drawnow;


% Generate threshold based off of AUC
[~,~,~,~,OPTROCPT] = perfcurve(labels,IIIlfS,posclass);
optTPR = OPTROCPT(1,2); %optimal TPR
[TPR,~,T,~] = perfcurve(labels,IIIlfS,posclass,'Xcrit','TPR');
ind = find(TPR==optTPR);
searchScore = mean(T(ind));
disp(strcat('Score for Patient 3:',num2str(searchScore)));

warning(s)  % restore the warning state


