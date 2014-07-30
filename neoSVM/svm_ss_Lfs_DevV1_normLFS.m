% SVM using searchScore and lowFreqScore
% requires clipScores

size(clipScores)

% your matlab & clips directory
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
nClips=length(clipDirNames);

% init vars for test subject
nthDir = 1; %choose test subject
clipLoc=strcat(matdir,clipsdir,clipDirNames{nthDir});
nIClips=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));
nIIClips=length(dir(strcat(clipLoc,filesep,'*_interictal_segment_*.mat')));
nIIIClips = nIClips + nIIClips;
lfS = clipScores{nthDir}(1,:);
sS = clipScores{nthDir}(2,:);

%% FOR DEBUGGING USING TEST SUBJECT

% Create true labels vector
labels = [ones(nIClips,1);zeros((nIIClips),1)];

%% Plot feature vector
figure;
scatter3(1:nIIIClips,lfS,sS);xlabel('Time');ylabel('LF Score'),zlabel('searchScore')
X=[lfS;sS]';

%% Train SVM
%svmStruct = svmtrain([lfS;sS],labels,'ShowPlot',true);
svmM = fitcsvm(X,labels,'ClassNames',[0 1],'Standardize',true);

cSvmM = compact(svmM);

%svmResult = svmclassify(svmStruct,[lfS;sS],'ShowPlot',true);
cSvmM = fitSVMPosterior(cSvmM,X,labels); %returns probabilities

% predict w/SVM
[predLabel,postProbs] = predict(cSvmM,X);
% figure;plot(predLabel);
figure;plot(postProbs(:,2))
qS = postProbs(:,2);

% normalize probabilities
% stdQ = std(qS(qS>0));
% uQ = mean(qS(qS>0));
% ceilQ = uQ+stdQ*2;
% qS(qS>ceilQ)= ceilQ; % set ceiling for outliers
% qS(qS>0)=(qS(qS>0)/ceilQ)/2+0.5; %normalize to ceiling b/n 0.5 to 1.0
% 
% stdQ = std(qS(qS<0));
% uQ = mean(qS(qS<0));
% ceilQ = uQ-stdQ*2;
% qS(qS<ceilQ)= ceilQ; % set ceiling for outliers
% qS(qS<0)=-(qS(qS<0)/ceilQ)/2; %normalize to ceiling b/n 0.0 to 0.5
% 
% qS(qS==0)=0.5;
% 
figure;plot(qS);title('SVM Posterior-Probabilities');set(gcf,'color','w');drawnow

% table(labels,predLabel,postProbs(:,2),'VariableNames',...
%    {'TrueLabels','PredictedLabels','PosteriorProbabilities'})

% Plot results
figure
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
    scatter3(i,lfS(i),sS(i),colr);xlabel('Time');ylabel('LF Score'),zlabel('searchScore');title('SVM results');hold on
end
drawnow;set(gcf,'color','w');

%% Compute ROC for SVM

% Tabulate scores
scores = qS;

% Define positive class
posclass = 1;

% Generate AUC
[X,Y,T,AUC] = perfcurve(labels,scores,posclass);

% Plot ROC Curve
figure;
plot(X,Y)
xlabel('False positive rate'); ylabel('True positive rate')
title('ROC, using mean landmarks per top 200 matches, basis X-axis')
annotation('textbox',...
    [0.3 0.65 0.2 0.1],...
    'String',{'AUC:',num2str(AUC)},...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','--',...
    'EdgeColor',[0 0 0],...
    'LineWidth',1,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);
axis square; set(gcf,'color','w');


break
%% Train SVM
% for real data, still uses same parameters
for ii=1:nClips
    lfS = clipScores{ii}(1,:);
    sS = clipScores{ii}(2,:);
    svmStruct = svmtrain([lfS,sS],labels);
    save(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/subject',num2str(ii)),svmStruct);
end

%%make submission
tic
resultsCell=cell(1,nClips);

for ii=1:nClips

    % load SVM classifier
    load(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/subject',num2str(ii)),svmStruct);    
    
    %without this line, hash tables from other clips might get used, comment out
    %for faster 1-file testing only
    clear get_hash_hits2
    clear find_landmarks2
    clear findLandmarks3
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    testClips=dir(strcat(clipLoc,filesep,'*_test_segment*.mat'));
    nTests=length(testClips);
    seizure=zeros(nTests,1);
    clip = cell(nTests,1);
    
     %ghetto progress
    disp(['finding submission values' clipLoc])
    
    for i=1:nTests
        
%         searchScore=sParLoad('mrscore.mat');
        R=match_query2(clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'));
        a=lfScoreTestSeg(clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'),ii);
        b=mean(R(:,2));
        
%         % flip 3 and 6 searchScores?
%         if strcmp(clipDirNames{ii},'Patient_3')||strcmp(clipDirNames{ii},'Patient_6')
% %             seizure(i,1)=mean(R(:,2))<searchScore(ii);
%             if a>b
%                 seizure(i,1)=1-(mean(R(:,2))/2/searchScore(ii)*(mean(R(:,2))<searchScore(ii)));
%             else
%                 seizure(i,1)=1-((1-searchScore(ii)/mean(R(:,2))/2)*(mean(R(:,2))>=searchScore(ii)));
%             end
%         else
% %             seizure(i,1)=mean(R(:,2))>searchScore(ii);
%             if a>b
%                 seizure(i,1)=mean(R(:,2))/2/searchScore(ii)*(mean(R(:,2))<searchScore(ii));
%             else
%                 seizure(i,1)=(1-searchScore(ii)/mean(R(:,2))/2)*(mean(R(:,2))>=searchScore(ii));
%             end
%         end

    % Classify w/SVM
    ScoreSVMModel = fitPosterior(SVMStruct,a,b); %returns probabilities
    seizure = ScoreSVMModel;
    
%         if SVM = seizure
%             %need to figure out how to get probability from SVM
%             seizure(i,1)=mean(R(:,2))/2/searchScore(ii)*(mean(R(:,2))<searchScore(ii));
%         else
%             %if not seizure
%             seizure(i,1)=(1-searchScore(ii)/mean(R(:,2))/2)*(mean(R(:,2))>=searchScore(ii));
%         end

        %create clip nmae
        clip{i}=strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat');
        
        
    end
    
    %assume early probability is same as seizure
    early=seizure;
    resultsCell{ii} = table(clip,seizure,early);
    save('sResults.mat','resultsCell');
    
end

submissionTable = vertcat(resultsCell{:});
writetable(submissionTable);

toc