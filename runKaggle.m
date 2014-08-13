%% Scipt for generating Kaggle submission
% load Tim's resultsCell
% this saves reuslts into SVMcell

global matdir
global clipsdir
global clipDirNames
global clipLoc
global scoMeth
global nIClips
global nIIClips
global nTSClips
global resultsCell

%your matlab directory & clips dir(relative to matlab, leave off^ that slash, as below. these get mushed together)
%matdir='C:\Users\paul\Documents\MATLAB\';
%clipsdir='kaggleShazam\clips\';
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';

%requires data to be in separate folders per patient with _ in folder name
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;

% best feature for each subject (prelim)
% fvOrder = {'mll','smmSco','rmsm','tmmSco','mll','mll','mspSco','rmsm','mll','rmsm','rmsm','rmsm'};

% Methods for each subject
%scoMeth{6} = {'mspSco 01,02,03,04','mspSco 13,14,15,16'}; % 'for P3 1-4 : 13-16
scoMeth{6} = {'mspSco 01,02,03,04','mspSco 13,14,15,16','mll 01,02,03,04','mll 13,14,15,16','rmsm 01,02,03,04','rmsm 13,14,15,16'}; % 'for P3 1-4 : 13-16
scoMeth{7} = {'mspSco0506','mspSco1011','rmsm 05,06','rmsm 10,11','phS'}; % 'for P3
scoMeth{10} = {'mspSco 07,08,09', 'mspSco 15,16', 'mspSco 23,24,25', 'rmsm'}; % 'for P6
    

tic
for ii=7:7 %set 1:12 to loop through all subjects
    %%load file data
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    nIClips=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));
    nIIClips=length(dir(strcat(clipLoc,filesep,'*_interictal_segment_*.mat')));
    nTSClips=length(dir(strcat(clipLoc,filesep,'*_test_segment_*.mat')));

    X = lfScoreTrainSVM(ii,ii);
    break
    %SVMlatencyTrain(X,i);
    
    size(resultsCell)

    %predsSVM = svmPredict(ii,resultsCell);
end

saveSubmission(predsSVM)

disp('Done! Total time spent:')
toc