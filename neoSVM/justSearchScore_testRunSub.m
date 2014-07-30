%% Loads previously generated testSegment searchScores into resultCell
% - use for mixing and matching seizure detection methods
% - P3: SVM 
% - P8: SVM

% your matlab & clips directory
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
nClips=length(clipDirNames);

%check existence of resultsCell
size(resultsCell)
%resultsSVMCell=cell(1,nClips);

tic

s = warning('off','all'); % turn all warnings off
if true,
    disp('Warnings off while running this code!') %lots of warnings from peakPoints not finding a peak
    warning('Here is a warning that doesn''t have an id.')
end

for ii=12:12 %use previously determined searchScores AUC of 0.79 
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
    if strcmp(clipDirNames{ii},'Patient_6')
        sS = 1-sS;
    end
    seizure = sS;
    
    for i=1:nTests    
        %create clip name
        clip{i}=strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat');
    end    
    
    %assume early probability is same as seizure
    early=seizure;
    predsSVM{ii} = table(clip,seizure,early);
    %save('test20140725_sSVMResults.mat','resultsSVMCell');
   toc 
end

