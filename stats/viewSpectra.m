

%your matlab directory & clips dir
% matdir='/Volumes/bobo/';
% clipsdir='kaggleClips/';
%matdir='C:\Users\paul\Documents\MATLAB\';

    matdir='/Users/idaniboy/Documents/MATLAB/';
    clipsdir='kaggleClips/ictalSegs/';

%clips dir(relative to matlab, leave off^ that slash, as below. these get mushed together)

%subject to run it on (do 1 at a time for now)
clipDirNames={'Patient_3'};

% %automation will require data to be in separate folders per patient with _ in folder name
% a=dir(strcat(matdir,clipsdir,'*_*'));
% clipLoc=strcat(matdir,clipsdir,clipDirNames{nthDir});
clipLoc=strcat(matdir,clipsdir,clipDirNames{1});

%type of segment to view
typeI='_ictal_segment';
% typeI='*_test_segment*';

%number of training segment randsamples to view 
nSamples=350;
a=dir(strcat(clipLoc,filesep,'*',typeI,'*.mat'));
srsOfSegments={a(randsample(length(a),min(nSamples,length(a)))).name};
nTests=length(srsOfSegments);

for i=1:nTests

    [data,freq,latency,channels]=sParLoad(strcat(clipLoc,filesep,srsOfSegments{i}));
    
    preLq= findLandmarksViewr(data,freq,channels);

end
