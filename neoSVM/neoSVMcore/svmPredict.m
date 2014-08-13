function predsSVM = svmPredict(subj,predsSVM)
%% Make submission from previously generated testSegment searchScores
% add p3 SVM

global matdir
global clipsdir
global clipDirNames
global clipLoc
global scoMeth
global nIClips
global nIIClips
global nTSClips

s = warning('off','all'); % turn all warnings off
if true,
    disp('Warnings off while running this code!') %lots of warnings from peakPoints not finding a peak
    warning('Here is a warning that doesn''t have an id.')
end

%% SVM

% load appropriate SVM classifier
load(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(subj)),'cSvmM');    
load(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Latency_Subject_',num2str(subj)),'cLatSvmM');    

TlfS = zeros(length(scoMeth{subj}),nTSClips);
clip = cell(nTSClips,1);

for nthMeth=1:length(scoMeth{subj});

    disp(['Finding ' scoMeth{subj}{nthMeth} ' searchScore for: ' clipLoc])

    % low Freq scores for ictal clips
    for i=1:nTSClips
        %create clip name
        clip{i,1}=strcat(clipDirNames{subj},'_test_segment_',num2str(i),'.mat');

        [data,~,freq,channels]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{subj},'_test_segment_',num2str(i),'.mat'));
        numCh = length(fieldnames(channels));

        % Get score
        TlfS(nthMeth,i) = lfScoreGetPerChan(data,scoMeth{subj}{nthMeth},numCh,freq,clipLoc);

    end
end
 
save(strcat('P',num2str(subj),'_TSscores.mat'),'TlfS')


%% post-processing
switch subj
    case 6
        % code for making ratios from first 3 feature pairs for P2
        Xtemp = TlfS; %X is longest vertically
        X = zeros(3,size(TlfS,2));
        X(1,:) = Xtemp(1,:)./Xtemp(2,:);
        X(2,:) = Xtemp(3,:)./Xtemp(4,:);
        X(3,:) = Xtemp(5,:)./Xtemp(6,:);
        X = X';
    otherwise
        X = TlfS';        
end

% join scores

%sS7 = 1-resultsCell{1, 7}.seizure(:,7);
%X = TlfS';

% predict w/SVM
[predLabel,postProbs] = predict(cSvmM,X);
seizure = postProbs(:,2);   

% predict seizureLatency <=15 sec w/SVM
%early=seizure; %assume early probability is same as seizure
[predLatLabel,postLatProbs] = predict(cLatSvmM,X);
early = postLatProbs(:,2);   

%save to table
predsSVM{subj} = table(clip,seizure,early);

save('sSVMResults20140806_try1.mat','predsSVM');
 
warning(s)  % restore the warning state
disp('done')




