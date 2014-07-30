%% Make submission from previously generated testSegment searchScores

% your matlab & clips directory
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
nClips=length(clipDirNames);

%check existence of resultsCell
size(resultsCell)
predsSVM=cell(1,nClips);
scoresSVM=cell(1,nClips);

tic

s = warning('off','all'); % turn all warnings off
if true,
    disp('Warnings off while running this code!') %lots of warnings from peakPoints not finding a peak
    warning('Here is a warning that doesn''t have an id.')
end


for ii=1:12 %use manualMode for patient 8
    tic
    % load appropriate SVM classifier
    load(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(ii)),'cSvmM');    
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    testClips=dir(strcat(clipLoc,filesep,'*_test_segment*.mat'));
    nTests=length(testClips);
    seizure=zeros(nTests,1);
    clip = cell(nTests,1);
    
     %ghetto progress
    disp(['finding submission values' clipLoc])

    %load old searchScores
    sS = cell2mat(table2cell(resultsCell{ii}(:,2)));
    sS = sS(:,ii);
    %flip pt 6
    if strcmp(clipDirNames{ii},'Patient_3')||strcmp(clipDirNames{ii},'Patient_6')
        sS = 1-sS;
    end
    
% 
%     lfS = zeros(nTests,1);
% 
%     % calculate lowFreqScores and make clipNames
%     bpD =  zeros(nTests,1);
%     bpB = zeros(nTests,1);
%     bpG = zeros(nTests,1);
    for i=1:nTests
        %lfS(i)=lfScoreTestSeg(ii,clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'));
        [bpD(i),bpB(i),bpG(i)]=bpScoreTestSeg(ii,clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'));
        %create clip name
        clip{i}=strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat');
    end
        
    figure
    plot(bpD,'r');hold on
    plot(bpB,'g');hold on
    plot(bpG,'b');hold on
    plot(lfS,'c');hold on
    plot(sS,'m');hold on
    

    % predict w/SVM
    [predLabel,postProbs] = predict(cSvmM,[bpD,bpB,bpG,sS]);
    seizure = postProbs(:,2);    

    
    %assume early probability is same as seizure
    early=seizure;
    predsSVM{ii} = table(clip,seizure,early);
    scoresSVM{ii} = [lfS,bpD,bpB,bpG];
    
    save('sSVMResults.mat','predsSVM','scoresSVM');
   toc 
end

%% 
for ii=12:12 %use manualMode for patient 8

    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    testClips=dir(strcat(clipLoc,filesep,'*_test_segment*.mat'));
    nTests=length(testClips);
    seizure=zeros(nTests,1);
    clip = cell(nTests,1);
    
    % init vars for test subject
    nIClips=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));
    nIIClips=length(dir(strcat(clipLoc,filesep,'*_interictal_segment_*.mat')));
    nIIIClips = nIClips + nIIClips;

    % Create true labels vector
    labels = [ones(nIClips,1);zeros((nIIClips),1)];    

    lfSTrain = clipScoresLFS{ii}(1,:);
    
    % Compute ROC for SVM
    posclass = 1;% Define positive class
    [X,Y,T,AUC] = perfcurve(labels,lfSTrain,posclass);% Generate AUC

%     % Plot ROC Curve
%     figure
%     plot(X,Y)
%     xlabel('False positive rate'); ylabel('True positive rate')
%     title(strcat('ROC, Subj_',num2str(ii)))
%     annotation('textbox',...
%         [0.9 0.1 0.1 0.1],...
%         'String',{'AUC:',num2str(AUC)},...
%         'FontSize',12,...
%         'FontName','Arial',...
%         'LineStyle','--',...
%         'EdgeColor',[0 0 0],...
%         'LineWidth',1,...
%         'BackgroundColor',[0.9  0.9 0.9],...
%         'Color',[0.84 0.16 0]);
%     axis square; set(gcf,'color','w');drawnow;    
%     save(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(ii)),'cSvmM');
%     drawnow;
    
    % Generate threshold based off of AUC
    [~,~,~,~,OPTROCPT] = perfcurve(labels,lfSTrain,posclass);
    optTPR = OPTROCPT(1,2); %optimal TPR
    [TPR,~,T,~] = perfcurve(labels,lfSTrain,posclass,'Xcrit','TPR');
    ind = find(TPR==optTPR);
    lfsThresh = mean(T(ind));
    disp(strcat('searchScore for Patient8:',num2str(lfsThresh)));

    % calculate lowFreqScores and make clipNames
    lfS = zeros(nTests,1);
    for i=1:nTests
        testSegLfs=lfScoreTestSeg(ii,clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'));
        if testSegLfs>=lfsThresh
            seizure(i,1)=1;
        else
            seizure(i,1)=(testSegLfs/lfsThresh)/2;
        end
        %create clip name
        clip{i}=strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat');
    end
        
    %assume early probability is same as seizure
    early=seizure;
    predsSVM{ii} = table(clip,seizure,early);
    
    save('sSVMResults.mat','predsSVM');
    
end

%% Use manual mode for patient 3


% your matlab & clips directory
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
nClips=length(clipDirNames);

% test subject
ii = 7;

% init vars for test subject
clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
nIClips=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));
nIIClips=length(dir(strcat(clipLoc,filesep,'*_interictal_segment_*.mat')));
nIIIClips = nIClips + nIIClips;
lfS = clipScoresP3{ii}(1,:);
%sS = clipScores{ii}(2,:);

testClips=dir(strcat(clipLoc,filesep,'*_test_segment*.mat'));
nTests=length(testClips);
seizure=zeros(nTests,1);
clip = cell(nTests,1);

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

% Generate threshold based off of AUC
[~,~,~,~,OPTROCPT] = perfcurve(labels,lfS,posclass);
optTPR = OPTROCPT(1,2); %optimal TPR
[TPR,~,T,~] = perfcurve(labels,lfS,posclass,'Xcrit','TPR');
ind = find(TPR==optTPR);
lfsThresh = mean(T(ind));
disp(strcat('searchScore for Patient3:',num2str(lfsThresh)));

s = warning('off','all'); % turn all warnings off
if true,
    disp('Warnings off while running this code!') %lots of warnings from peakPoints not finding a peak
    warning('Here is a warning that doesn''t have an id.')
end

    % calculate lowFreqScores and make clipNames
    lfS_SS = zeros(nTests,1);
    lfS_phS = zeros(nTests,1);

    for i=1:nTests
        %lfS_SS(i)=lfScoreTestSeg2(ii,clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'));
        lfS_phS(i)=lfScoreTestSeg2(ii,clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'));

        %create clip name
        clip{i}=strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat');
    end
    
    %% stats
    % load slopeScoresp3... 
    figure;plot(lfS_SS); hold on
    plot(1,lfsThresh,nTests,lfsThresh,'r');
    
    figure;plot(lfS)
    
    %% Build table
    lfsThresh = 1.1;
    maxSS = max(lfS_SS);
    ceilSS = uSS+2*stdSS;
    stdSS = std(lfS_SS);
    uSS = mean(lfS_SS);
    norm_lfS_SS = lfsThresh;
    norm_lfS_SS(norm_lfS_SS>ceilSS) = ceilSS;
    for i=1:nTests
        if lfS_SS(i)>=lfsThresh
            seizure(i,1)=((maxSS - lfS_SS(i)-lfsThresh)/(maxSS - lfsThresh))/2+0.49;
        else
            seizure(i,1)=(lfS_SS(i)/lfsThresh)/2;
        end
    end
    warning(s)  % restore the warning state
    
    %assume early probability is same as seizure
    early=seizure;
    resultsSVMCell{ii} = table(clip,seizure,early);
    
    save('sSVMResults.mat','predsSVM');
    disp('done');