function createHashTable()
%best to do one at a time for now, fill in info below

    %your matlab directory
    matdir='/Users/idaniboy/Documents/MATLAB/';
    %matdir='C:\Users\paul\Documents\MATLAB\';
    %clips dir(relative to matlab, leave off^ that slash, as below. these get mushed together)
    %clipsdir='kaggleShazam\clips\';
    clipsdir='kaggleClips/ictalSegs/';

    %requires data to be in separate folders per patient with _ in folder name
    a=dir(strcat(matdir,clipsdir,'*_*'));
    clipDirNames={'Dog_1'};%{a.name};%{'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
    
    for nthDir=1:length(clipDirNames);
        clear save_hashes
        clipLoc=strcat(matdir,clipsdir,clipDirNames{nthDir});
        %load concatenated_ filenames
        nSeizures=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));
%         %used to track which seizure within a folder is currently active
%         aSHIndex=1;
        H=[];
        for nthSeizure=1:nSeizures
            [data,freq,~,channels]=sParLoad(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_segment_',num2str(nthSeizure),'.mat'));
            disp('total # of seizures:')
            disp(nSeizures)
            %preallocate and group all seizures in same patient folder seperately
            preH = landmark2hash(findLandmarks0(data,ceil(freq),channels));
            
            %might work speeding up with parfor, needs work
%             nLandmarks=size(H,1);
%             if nthSeizure==1
%                 allSeizureHashes=[H;zeros(round(nSeizures*1.3*nLandmarks),3)];
%                 aSHIndex=aSHIndex+nLandmarks;
%             elseif nthSeizure==nSeizures
%                 allSeizureHashes(aSHIndex:aSHIndex+nLandmarks-1,:)=H;
%                 allSeizureHashes=allSeizureHashes(1:aSHIndex+nLandmarks-1,:);
%             else
%                 allSeizureHashes(aSHIndex:aSHIndex+nLandmarks-1,:)=H;
%             end
            H=[H; preH];
%             nthSeizure
        end
        save_hashes(H,clipLoc,clipDirNames{nthDir});
        nthDir
        %for parfor if ever necessary
%         save_hashes(allSeizureHashes,clipLoc,clipDirNames{nthDir});
        
    end

end


