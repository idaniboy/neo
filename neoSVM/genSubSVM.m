%% Make submission with SVM
% derived from svm_ss_Lfs_TestRun
% calculates searchScore and lowFreqScore from testSegs

tic
resultsCell=cell(1,nClips);

for ii=1:11 %use manualMode for patient 8

    % load SVM classifier
    load(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(ii)),cSvmM);    
    
    %without this line, hash tables from other clips might get used, comment out
    %for faster 1-file testing only
    clear get_hash_hits2
    clear find_landmarks2
    clear findLandmarks3
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    testClips=dir(strcat(clipLoc,filesep,'*_test_segment*.mat'));
    nTests=length(testClips);
    seizure=zeros(nTests,1);
    clip = cell(nTests,1);
    
     %ghetto progress
    disp(['finding submission values' clipLoc])
    
    for i=1:nTests
        
%         searchScore=sParLoad('mrscore.mat');
        R=match_query2(clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'));
        lfS=lfScoreTestSeg(clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'),ii);
        b=mean(R(:,2));
        
        %flip pt 3 and 6
        if strcmp(clipDirNames{ii},'Patient_3')||strcmp(clipDirNames{ii},'Patient_6')
%             seizure(i,1)=mean(R(:,2))<searchScore(ii);
            if a>b
                sS=1-(mean(R(:,2))/2/searchScore(ii)*(mean(R(:,2))<searchScore(ii)));
            else
                sS=1-((1-searchScore(ii)/mean(R(:,2))/2)*(mean(R(:,2))>=searchScore(ii)));
            end
        else
%             seizure(i,1)=mean(R(:,2))>searchScore(ii);
            if a>b
                sS=mean(R(:,2))/2/searchScore(ii)*(mean(R(:,2))<searchScore(ii));
            else
                sS=(1-searchScore(ii)/mean(R(:,2))/2)*(mean(R(:,2))>=searchScore(ii));
            end
        end
        
        % predict w/SVM
        [predLabel,postProbs] = predict(cSvmM,[lfS;sS]');
        seizure = postProbs(:,2);    
    
        %create clip name
        clip{i}=strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat');
        
    end
    
    %assume early probability is same as seizure
    early=seizure;
    resultsCell{ii} = table(clip,seizure,early);
    save('sResults.mat','resultsCell');
    
end

submissionTable = vertcat(resultsCell{:});
writetable(submissionTable);

toc