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
            
            if TSlfS(i,nthMeth)>2
                TSlfS(i,nthMeth)=2;
            end

            if TSlfS(i,nthMeth)>=1
                seizure = TSlfS(i,nthMeth)/2;
            else
                seizure = ((TSlfS(i,nthMeth)-0.4142)/(1-0.4142))/2;
            end
        end
        
    end
end



%% assume early probability is same as seizure
early=seizure;
predsSVM{ii} = table(clip,seizure,early);
scoresSVM{ii} = X;
