%% Make submission from previously generated testSegment searchScores
% add p3 SVM

% your matlab & clips directory
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
nClips=length(clipDirNames);

%check existence of resultsCell
size(resultsCell)
%predsSVM=cell(1,nClips);
%scoresSVM=cell(1,nClips);

tic

s = warning('off','all'); % turn all warnings off
if true,
    disp('Warnings off while running this code!') %lots of warnings from peakPoints not finding a peak
    warning('Here is a warning that doesn''t have an id.')
end

%% use SVM for patient 3
for ii=7:7 
    
    % load appropriate SVM classifier
    load(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(ii)),'cSvmM');    
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    testClips=dir(strcat(clipLoc,filesep,'*_test_segment*.mat'));
    nTests=length(testClips);
    scoMeth = {'mspSco0506','mspSco1011','rmsmT3','phS'};

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
             TlfS(nthMeth,i) = lfScoreGet(data,scoMeth{nthMeth},numCh,freq,clipLoc);

        end
    end
end   

%% join scores

%sS7 = 1-resultsCell{1, 7}.seizure(:,7);
%X = TlfS';
X = [TlfS(1,:);TlfS(2,:);TlfS(3,:)]';

% predict w/SVM
[predLabel,postProbs] = predict(cSvmM,X);
seizure = postProbs(:,2);   
    
%assume early probability is same as seizure
early=seizure;
predsSVM{ii} = table(clip,seizure,early);
scoresSVM{ii} = X;

save('sSVMResults20140802_try1.mat','predsSVM','scoresSVM');
toc 
warning(s)  % restore the warning state
disp('done')




