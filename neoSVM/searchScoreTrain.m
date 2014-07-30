% % creates feature vectors from high frequency features of seizure data
% save searchScores for each clip on training data to train SVM
% run lfScore first to make clipScores array!
% this version was made to use mean searchScores to determine seizure prob

%your matlab directory
if ~ispc
    % your matlab & clips directory
    matdir='/Volumes/bobo/';
    clipsdir='kaggleClips/';
%     matdir='/Users/idaniboy/Documents/MATLAB/';
%     clipsdir='kaggleClips/ictalSegs/';
    a=dir(strcat(matdir,clipsdir,'*_*'));
a(1) = [];  %remove 1st element in mac because it is .ds_store
    
else
    matdir='C:\Users\paul\Documents\MATLAB\';
    clipsdir='kaggleShazam\clips\';
    % %automation requires data to be in separate folders per patient with _ in folder name
    a=dir(strcat(matdir,clipsdir,'*_*'));
end

clipDirNames={a.name};
% clipDirNames={'Patient_1'};
% clipDirNames={'Dog_1'};%,'Dog_2','Dog_3'};%,'Dog_2','Dog_3'};
nClips=length(clipDirNames);

% load searchScore values
searchScore=sParLoad('mrscore.mat');


%check clipScores
disp(strcat('Size clipScores:',num2str(size(clipScores))));

%find searchScore
tic
for ii=1:nClips
    
    %without this line, hash tables from other clips might get used, comment out
    %for faster 1-subject testing only
    clear get_hash_hits2
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});        
    nIClips=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));    
    nIIClips=length(dir(strcat(clipLoc,filesep,'*_interictal_segment_*.mat')));
    
    totalI=zeros(1,nIClips);
    totalII=zeros(1,nIIClips);
    
    a=searchScore(ii);    
    
    %ghetto progress 
    disp(['finding searchScore for ' clipLoc])
   

    
    %roundabout way to allow for mean/std in parfor
    parfor i=1:nIClips
        
        R=match_query2(clipLoc,strcat(clipDirNames{ii},'_ictal_segment_',num2str(i),'.mat'));
                b=mean(R(:,2));

%         if ~isempty(R)
%             totalI(i)=R(1,2);
%         end
%             total(i)=mean(R(:,2));

        if strcmp(clipDirNames{ii},'Patient_3')||strcmp(clipDirNames{ii},'Patient_6')
%             seizure(i,1)=mean(R(:,2))<searchScore(ii);
            if a>b
                totalI(i)=1-(mean(R(:,2))/2/searchScore(ii)*(mean(R(:,2))<searchScore(ii)));
            else
                totalI(i)=1-((1-searchScore(ii)/mean(R(:,2))/2)*(mean(R(:,2))>=searchScore(ii)));
            end
        else
%             seizure(i,1)=mean(R(:,2))>searchScore(ii);
            if a>b
                totalI(i)=mean(R(:,2))/2/searchScore(ii)*(mean(R(:,2))<searchScore(ii));
            else
                totalI(i)=(1-searchScore(ii)/mean(R(:,2))/2)*(mean(R(:,2))>=searchScore(ii));
            end
        end

    end
    
    parfor i=1:nIIClips
        R=match_query2(clipLoc,strcat(clipDirNames{ii},'_interictal_segment_',num2str(i),'.mat'));
                b=mean(R(:,2));

%         if ~isempty(R)
%             totalII(i)=R(1,2);
%         end

        if strcmp(clipDirNames{ii},'Patient_3')||strcmp(clipDirNames{ii},'Patient_6')
%             seizure(i,1)=mean(R(:,2))<searchScore(ii);
            if a>b
                totalII(i)=1-(mean(R(:,2))/2/searchScore(ii)*(mean(R(:,2))<searchScore(ii)));
            else
                totalII(i)=1-((1-searchScore(ii)/mean(R(:,2))/2)*(mean(R(:,2))>=searchScore(ii)));
            end
        else
%             seizure(i,1)=mean(R(:,2))>searchScore(ii);
            if a>b
                totalII(i)=mean(R(:,2))/2/searchScore(ii)*(mean(R(:,2))<searchScore(ii));
            else
                totalII(i)=(1-searchScore(ii)/mean(R(:,2))/2)*(mean(R(:,2))>=searchScore(ii));
            end
        end
        
    end
    
    %needed for some reason, otherwise the 2nd parfor loop considers total a temporary variable instead of sliced
    %    total(ii,:,2)=bstemp;
    
    
    
    clipScores{ii}(2,:) = [totalI,totalII]; %create placeholder for searchScore
    
end

toc

save('cScores.mat','clipScores')
