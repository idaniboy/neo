% function createHashTable
% derived from bpSD2forClips
% creates vector of delta, alpha/beta, and gamma band powers
% twice as fast as bpSDforClips.m on Dog1
% AUC is ~ the same
% also generates submission file (currently broken but tim fixing)

    %your matlab directory
%     matdir='/Users/idaniboy/Documents/MATLAB/';
%     clipsdir='kaggleClips/ictalSegs/';
    matdir='/Volumes/bobo/';
    clipsdir='kaggleClips/';

    %requires data to be in separate folders per patient with _ in folder name
    a=dir(strcat(matdir,clipsdir,'*_*'));
    %clipDirNames={'Dog_1'};%{a.name};
    clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;

    % clipDirNames={'Dog_2'};%,'Dog_2','Dog_3'};
    nClips=length(clipDirNames);

%% find bpScore
tic
for nthDir=1:1;

    clipLoc=strcat(matdir,clipsdir,clipDirNames{nthDir});

    %load concatenated_ filenames
    nIClips=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));
    nIIClips=length(dir(strcat(clipLoc,filesep,'*_interictal_segment_*.mat')));
    nTSClips=length(dir(strcat(clipLoc,filesep,'*_test_segment_*.mat')));

    %initialize some parameters
    [~,~,freq,channels]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_segment_',num2str(1),'.mat'));
    freq = ceil(freq);
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
    end
    numCh = length(fieldnames(channels));
    data = zeros(numCh,(size(zB,2)*2+freq));

    % initialize mean bandpower per 1sec clip
    mpID = zeros(1,nIClips);
    mpIB = zeros(1,nIClips);
    mpIG = zeros(1,nIClips);
    mpIV = zeros(1,nIClips);
    
    disp(['Finding bpS for: ' clipLoc])

    % bandpower for ictal clips
    for i = 1:nIClips;
        [datum,~,~,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_segment_',num2str(i),'.mat'));

        switch freq
            case 5000
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
        ppV = zeros(numCh,1);
         
        for j = 1:size(data,1)
            %Bandpower for selected band at each channel, j
            ppD(j) = bandpower(data(j,:),freq,[1 10]);
            ppB(j) = bandpower(data(j,:),freq,[20 55]);
            ppG(j) = bandpower(data(j,:),freq,[40 topFreq]);
            
            %Variance
            %ppV(j) = var(data(j,:),[],2);
        end

        % Mean bandpower at each 1sec clip
        mpID(i) = mean(ppD);
        mpIB(i) = mean(ppB);        
        mpIG(i) = mean(ppG);
        mpIV(i) = mean(ppV);
    end

    mpIID = zeros(1,nIIClips);
    mpIIB = zeros(1,nIIClips);
    mpIIG = zeros(1,nIIClips);
    mpIIV = zeros(1,nIIClips);

    % bandpower for interictal clips
    for i = 1:nIIClips;
        [datum,~,~,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_interictal_segment_',num2str(i),'.mat'));
        switch freq
            case 5000
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
        ppV = zeros(numCh,1);
        
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
        mpIIV(i) = mean(ppV);
    end

    % combine ictal and interictal
    mpD = [mpID,mpIID];
    mpB = [mpIB,mpIIB];
    mpG = [mpIG,mpIIG];
    mpV = [mpIV,mpIIV];

    % entire segment statistics
    UG = mean(mpG);     
    stdUG = std(mpG); %std dev
    medAllUG = median(mpG); %median
    
    
    %% Train SVM

    lfS = clipScoresLFS{1}(1,:);
    sS = clipScores{1}(2,:);    

     %% Plot
    figure
    plot(mpD,'b');hold on;
    plot(mpB,'g');hold on;
    plot(mpG,'r');hold on;   
    plot(mpV,'c');hold on;
    
    % Training data
    X=[mpD;mpB;mpG;sS]';

    % Create true labels vector
    labels = [ones(nIClips,1);zeros((nIIClips),1)];    
    
    svmM = fitcsvm(X,labels,'ClassNames',[0 1],'Standardize',true,'Cost',[0,1;2,0]); 
    cSvmM = compact(svmM);
    cSvmM = fitSVMPosterior(cSvmM,X,labels); %returns probabilities

    % predict w/SVM
    [predLabel,postProbs] = predict(cSvmM,X);
    
    % plot Posterior Probabilities
    qS = postProbs(:,2);    
    figure;subplot(2,1,1)
    plot(qS);title('SVM Posterior-Probabilities');set(gcf,'color','w');drawnow
    
    % Plot prediction results
%     figure
%     for i=1:size(predLabel,1)
%         if predLabel(i) ~= labels(i)
%             colr = 'r';
%         else
%             if labels(i) == 1
%                 colr = 'g';
%             else
%                 colr = 'b';
%             end
%         end
%         scatter3(i,lfS(i),sS(i),colr);xlabel('Time');ylabel('LF Score'),zlabel('searchScore');title('SVM results');hold on
%     end
%     drawnow;set(gcf,'color','w');
    
    %% Compute ROC 
    % Tabulate scores
    scores = qS;

    % Define positive class
    posclass = 1;

    % Generate AUC
    [X,Y,T,AUC] = perfcurve(labels,scores,posclass);

    % Plot ROC Curve
    subplot(2,1,2)
    plot(X,Y)
    xlabel('False positive rate'); ylabel('True positive rate')
    title('ROC, using mean landmarks per top 200 matches, basis X-axis')
    annotation('textbox',...
        [0.9 0.1 0.2 0.1],...
        'String',{'AUC:',num2str(AUC)},...
        'FontSize',14,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[0 0 0],...
        'LineWidth',1,...
        'BackgroundColor',[0.9  0.9 0.9],...
        'Color',[0.84 0.16 0]);
    axis square; set(gcf,'color','w');

    % Generate gamma score
    [~,~,~,~,OPTROCPT] = perfcurve(labels,scores,posclass);
    [TPR,~,T,~] = perfcurve(labels,scores,posclass,'Xcrit','TPR');
    gammaScore(nthDir) = mean(T(TPR==OPTROCPT(1,2)));
    disp(strcat('Gamma score: ',num2str(gammaScore(nthDir))));
    toc

end



