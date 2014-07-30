%% SVM_insertP8

% SVM
for ii=12:12 
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
    %flip pt 3 and 6
    if strcmp(clipDirNames{ii},'Patient_3')||strcmp(clipDirNames{ii},'Patient_6')
        sS = 1-sS;
    end
    

    lfS = zeros(nTests,1);

    % calculate lowFreqScores and make clipNames
    bpD =  zeros(nTests,1);
    bpB = zeros(nTests,1);
    bpG = zeros(nTests,1);
    for i=1:nTests
        lfS(i)=lfScoreTestSeg(ii,clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'));
        [bpD(i),bpB(i),bpG(i)]=bpScoreTestSeg(ii,clipLoc,strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat'));
        %create clip name
        clip{i}=strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat');
    end
        
    figure
    plot(bpD,'r');hold on
    plot(bpB,'g');hold on
    plot(bpG,'b');hold on
%    plot(lfS,'c');hold on
    plot(sS,'m');hold on
    

    % predict w/SVM
    [predLabel,postProbs] = predict(cSvmM,[bpD,bpB,bpG,sS]);
    seizure = postProbs(:,2);    

    
    %assume early probability is same as seizure
    early=seizure;
    resultsSVMCell{ii} = table(clip,seizure,early);
    scoresSVM{ii} = [lfS,bpD,bpB,bpG];
    
    %save('sSVMResults.mat','predsSVM');
   toc 
end

