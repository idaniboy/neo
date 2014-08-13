%% Make submission from previously generated testSegment searchScores
% add p6 SVM
% also includes latency SVM

% your matlab & clips directory
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
nClips=length(clipDirNames);

%check existence of resultsCell
%size(resultsCell)
%predsSVM=cell(1,nClips);
%scoresSVM=cell(1,nClips);

tic

s = warning('off','all'); % turn all warnings off
if true,
    disp('Warnings off while running this code!') %lots of warnings from peakPoints not finding a peak
    warning('Here is a warning that doesn''t have an id.')
end

%% use SVM for patient 3
for ii=10:10 
    
    % load appropriate SVM classifier
    load(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(ii)),'cSvmM');    
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    testClips=dir(strcat(clipLoc,filesep,'*_test_segment*.mat'));
    nTests=length(testClips);
    scoMeth = {'mspSco 07,08,09', 'mspSco 15,16', 'mspSco 23,24,25', 'rmsm'}; % 'for P6

    %data stores
    TlfS = zeros(3,nTests); %ictal low freq score    
    seizure=zeros(nTests,1);
    clip = cell(nTests,1);
    
    for nthMeth=1:length(scoMeth);
        
        disp(['Finding ' scoMeth{nthMeth} ' searchScore for: ' clipLoc])

        % low Freq scores for ictal clips
        for i=1:nTests
            %create clip name
            clip{i,1}=strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat');

            [data,~,freq,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_test_segment_',num2str(i),'.mat'));

             % Get score
             TlfS(nthMeth,i) = lfScoreGetPerChan(data,scoMeth{nthMeth},numCh,freq,clipLoc);

        end
    end
end   

%% join scores
load(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(ii)),'cSvmM');
%sS10 = 1-resultsCell{1, 10}.seizure(:,10);
%X = TlfS';
X = [TlfS(1:3,:)]';

% predict w/SVM
[predLabel,postProbs] = predict(cSvmM,X);
seizure = postProbs(:,2);   
disp('done predicting')
figure;plot(seizure)

% predict seizureLatency <=15 sec w/SVM
%early=seizure; %assume early probability is same as seizure
[predLatLabel,postLatProbs] = predict(cLatSvmM,X);
early = postLatProbs(:,2);   

%store table into predsSVM
predsSVM{ii} = table(clip,seizure,early);
scoresSVM{ii} = X;

save('sSVMResults20140803.mat','predsSVM','scoresSVM');
toc 
warning(s)  % restore the warning state
disp('done')




