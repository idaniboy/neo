% % creates feature vectors from low frequency features of seizure data
% development for kaggle seizure contest...
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

% initialize clipScoresrray to score lfScores
clipBpScores = cell(12,1);
    
    
%% find BpScore,
tic

for nthDir=1:nClips;  
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{nthDir});

    %load concatenated_ filenames
    nIClips=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));
    nIIClips=length(dir(strcat(clipLoc,filesep,'*_interictal_segment_*.mat')));
    nTSClips=length(dir(strcat(clipLoc,filesep,'*_test_segment_*.mat')));

    mpIG = zeros(1,nIClips);

    %initialize some parameters
    [~,freq,~,channels]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_segment_',num2str(1),'.mat'));
    freq = ceil(freq);
    numCh = length(fieldnames(channels));

    switch freq
        case 400
            zB = zeros(1,312); % total length = 1024
            topFreq = 200;
        case 500
            zB = zeros(1,262); % total length = 1024
            topFreq = 250;
        case 5000
            zb = zeros(1,12);% total length = 1024
            topFreq = 400;
            freq = 1000;
    end
    numCh = length(fieldnames(channels));
    data = zeros(numCh,(size(zB,2)*2+freq));

    % initialize mean bandpower per 1sec clip
    mpID = zeros(1,nIClips);
    mpIB = zeros(1,nIClips);
    mpIG = zeros(1,nIClips);
    mpIV = zeros(1,nIClips);
    
    disp(['Finding bpS for: ' clipLoc])

    % low Freq scores for ictal clips
    for i = 1:nIClips;
        [datum,~,~,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_segment_',num2str(i),'.mat'));
        switch freq
            case 1000
                %downsample freq from 5000 --> 1000;
                dfreq = 5;
                Ds = zeros(size(datum,1),size(datum,2)/5); %assumes SR 5000
                for j=1:size(datum,1)
                    Ds(j,:) = decimate(datum(j,:),dfreq);
                end
                datum = Ds;
        end

        data(:,zB+1:zB+freq) = datum;
        
        ppD = zeros(numCh,1);
        ppB = zeros(numCh,1);
        ppG = zeros(numCh,1);
        %ppV = zeros(numCh,1);
        
        for j = 1:numCh

            %Bandpower for selected band at each channel, j
            ppD(j) = bandpower(data(j,:),freq,[1 10]);
            ppB(j) = bandpower(data(j,:),freq,[20 55]);
            ppG(j) = bandpower(data(j,:),freq,[40 topFreq]);            
           
        end

        % Mean bandpower at each 1sec clip
        mpID(i) = mean(ppD);
        mpIB(i) = mean(ppB);        
        mpIG(i) = mean(ppG);
        %mpIV(i) = mean(ppV);
        
    end
    
   
    mpIID = zeros(1,nIIClips);
    mpIIB = zeros(1,nIIClips);
    mpIIG = zeros(1,nIIClips);
    %mpIIV = zeros(1,nIIClips);

    % bandpower for interictal clips
    for i = 1:nIIClips;
        [datum,~,~,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_interictal_segment_',num2str(i),'.mat'));
        switch freq
            case 1000
                %downsample if freq is 5000
                dfreq = 5;
                Ds = zeros(size(datum,1),size(datum,2)/5); %assumes SR 5000
                for j=1:size(datum,1)
                    Ds(j,:) = decimate(datum(j,:),dfreq);
                end
                datum = Ds;
        end
        data(:,zB+1:zB+freq) = datum;
        
        ppD = zeros(numCh,1);
        ppB = zeros(numCh,1);
        ppG = zeros(numCh,1);
        %ppV = zeros(numCh,1);
        
        for j = 1:size(data,1)
            %Bandpower for selected band at each channel, j
            ppD(j) = bandpower(data(j,:),freq,[1 10]);
            ppB(j) = bandpower(data(j,:),freq,[20 55]);
            ppG(j) = bandpower(data(j,:),freq,[40 topFreq]);
            %Variance
            %ppV(j) = var(data(j,:),[],2);
        end

        % Mean bandpower at each 1sec clip
        mpIID(i) = mean(ppD);
        mpIIB(i) = mean(ppB);   
        mpIIG(i) = mean(ppG);
        %mpIIV(i) = mean(ppV);
    end

    % combine ictal and interictal
    mpD = [mpID,mpIID];
    mpB = [mpIB,mpIIB];
    mpG = [mpIG,mpIIG];
    %mpV = [mpIV,mpIIV];

    %normalize and save    
    clipBpScores{nthDir}(1,:) = mpD; % low freq score in row 1 of clipScores
    clipBpScores{nthDir}(2,:) = mpB; % normalized low freq score in row 1 of clipScores
    clipBpScores{nthDir}(3,:) = mpG; %create placeholder for searchScore
    clipBpScores{nthDir}(4,:) = zeros(1,size(clipScores{nthDir},2)); %create placeholder for searchScore
    
    
end



save('cBpScores.mat','clipBpScores')


