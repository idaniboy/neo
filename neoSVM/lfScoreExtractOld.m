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
        scoMeth = {'mspSco','phS'};
        
        %data store for lfS
        lfS = zeros(nIClips+nIIClips,length(scoMeth));
        
    for nthMeth=1:length(scoMeth);
        
        %initialize prescores
        mspSco = zeros(numCh,1); % 
        sspSco = zeros(numCh,1);
        aspSco = zeros(numCh,1);
        mmpSco = zeros(numCh,1);
        tmpSco = zeros(numCh,1);
        rmsSco = zeros(numCh,1); %RMS pre-score    
        rmsT3Sco = zeros(numCh,1); %RMS pre-score    
        varS = zeros(numCh,1);
        phS = zeros(numCh,1); %phS = phase score
        rmsIZ = zeros(numCh,1); %phS = ictalZone RMS

        %data stores
        IlfS = zeros(1,nIClips); %ictal low freq score
        IIlfS = zeros(1,nIIClips); %interictal low freq score

        disp(['Finding searchScore for: ' clipLoc])

        % low Freq scores for ictal clips
        for i = 1:nIClips;
            [data,~,freq,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_segment_',num2str(i),'.mat'));
             switch freq
                case 5000
                    %downsample if freq is 5000
                    dfreq = 5;
                    Ds = zeros(size(data,1),size(data,2)/5); %assumes SR 5000
                    for j=1:size(data,1)
                        Ds(j,:) = decimate(data(j,:),dfreq);
                    end
                    data = Ds;
             end
             
             %Pre-process
             switch scoMeth{nthMeth} 
                case 'phS'
                        data = medfilt1(data',20)';
                otherwise
                    for j = 1:numCh
                        data(j,:) = lowPassFilterTest(data(j,:),freq/dfreq);
                    end   
             end

            for j = 1:numCh
                %select feature for feature vectors
                switch scoMeth{nthMeth} %subject number = nthDir

                    case 'rmsm'
                    % calculate RMS
                    rmsP = rms(data(j,:), 2); 
                    rmsSco(j) = rmsP;
                    
                    case 'rmsmT3'
                    % calculate RMS average for Top3 / RMS average for all
                    rmsT3P = rms(data(j,:), 2); 
                    if j == 5 || j == 6
                        rmsSco(j) = rmsT3P;
                    else
                        rmsSco(j) = rmsT3P;
                    end
                    rmsT3Sco(j) = rmsT3P;     

                    case 'mll'
                        % calculate lineLength
                        ll(j) = lineLength(data(j,:));

                    case 'var'
                       varS(j)=var(data(j,:));

                    case {'tmmSco','smmSco','mspSco'}

                        %peakFinding
                        [pks,locs] = findpeaks(data(j,:),'minpeakdistance',100,'MINPEAKHEIGHT',5);
                        if numel(locs) ~= 0
                            mLocs = max(locs);
                        else
                            mLocs = 0.1;
                        end
                        [npks,nlocs] = findpeaks(mLocs-data(j,:),'minpeakdistance',100,'MINPEAKHEIGHT',5);
                        npks = mLocs-npks;
                        % %plot
    %                    figure;plot(data(j,:)), hold on; plot(locs,pks+0.05,'k^','markerfacecolor',[0 1 1]); hold on; plot(nlocs,npks+0.05,'k^','markerfacecolor',[1 0 0]), hold off
                        %check if peaks founds
                        if numel(nlocs) == 0 
                            switch scoMeth{nthMeth}
                                case 'smmSco'
                                    % no peaks so zero score
                                    mmpSco(j) = 0;
                                case 'tmmSco'
                                    % no peaks so zero score
                                    tmpSco(j) = 0;
                                case 'mspSco'
                                    mspSco(j) = 0;
                                    %sspSco(j) = 0;
                            end
                        else

                            switch scoMeth{nthMeth}
                                case 'smmSco'
                                    % calculate sum max/min score
                                    mmpSco(j) = size(pks,2);
                                case 'tmmSco'
                                    % calculate thresholded max/min score
                                    tPks = pks(pks>100);
                                    tmpSco(j) = size(tPks,2);
                                case 'mspSco'
                         
                                    % calculate slopeScore, sS, parameters
                                    a = [pks',locs'];
                                    b = [npks',nlocs'];
                                    sS = vertcat(a,b);
                                    sS = sortrows(sS,2);

                                    % check if slope calcuable 
                                    if numel(sS) > 1
                                        % calculate slopeScore, sS
                                        sS2 = circshift(sS,-1); 
                                        sSS = abs((sS2(:,1) - sS(:,1)) ./ (sS2(:,2) - sS(:,2)));% calc dY / dX
                                        mspSco(j) = mean(sSS);% calculate mean slope score
                                        %sspSco(j) = sum(sSS);% calculate sum slope score
                                    else
                                        mspSco(j) = 0;
                                        %sspSco(j) = 0;
                                    end
                            end
                        end

                    case 'mAs'
                        % calculate mean amplitude score
                        n = 100; %number of max values to find
                        [aa,ix] = sort(data(j,:),'descend');
                        aa = aa(1:n);
                        % ix = ix(1:n);
                        aspSco(j) = sum(aa); 
                end

            end

            % Mean score at each 1sec clip
            switch scoMeth{nthMeth} %subject number = nthDir
                    case 'rmsm'
                        IlfS(i) = mean(rmsSco); 
                    case 'rmsmT3'
                        %find top 7 values
                        A = rmsT3Sco;
                        Asorted = sort(A(:),'descend');
                        IlfS(i)=mean(Asorted(1:7))/mean(Asorted(8:end));
                    case 'rmsIZ'           
                        rms5 = rms(data(5,:), 2); 
                        rms6 = rms(data(6,:), 2); 
                        rms9 = rms(data(6,:), 2); 
                        rms10 = rms(data(6,:), 2); 

                        IlfS(i)=rms5+rms6+rms9+rms10;                     
                    case 'phS'
                        IlfS(i)=sum((data(5,:))-data(6,:));
                         %IlfS(i)=sum(abs(data(5,:))+abs(data(6,:)));

                    case 'mll'
                        IlfS(i) = mean(ll);
                    case 'var'
                        %find top 5 values
    %                     A = varS;
    %                     Asorted = sort(A(:),'descend');
    %                     IlfS(i)=mean(Asorted(1:5))/mean(Asorted(5:end));                   
                          IlfS(i)=varS(5)+varS(6)+varS(9)+varS(10);                   

                    case 'smmSco'

                        A = mmpSco;
                        Asorted = sort(A(:),'descend');
                        IlfS(i)=mean(Asorted(1:7))/mean(Asorted(8:end));                    
                        %IlfS(i) = mean(mmpSco); % max/min score
                    case 'tmmSco'           
                        A = tmpSco;
                        Asorted = sort(A(:),'descend');
                        IlfS(i)=mean(Asorted(1:7))/mean(Asorted(8:end)); 
                        
                    case 'mspSco'
                        IlfS(i) = mean(mspSco);

            end

    %         IlfS(i) = mean(sspSco);
    %         IlfS(i) = mean(aspSco);


        end

        % mpIIMG = zeros(1,nIIClips);
        % mpIIHG = zeros(1,nIIClips);

        %initialize prescores
        mspSco = zeros(numCh,1); % 
        sspSco = zeros(numCh,1);
        aspSco = zeros(numCh,1);
        mmpSco = zeros(numCh,1);
        tmpSco = zeros(numCh,1);
        rmsSco = zeros(numCh,1); %RMS pre-score  
        rmsT3Sco = zeros(numCh,1);
        varS = zeros(numCh,1);

        % bandpower for interictal clips
        for i = 1:nIIClips;
            [data,~,~,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_interictal_segment_',num2str(i),'.mat'));

            %Pre-process
             switch scoMeth{nthMeth} 
                case 'phS'
                        data = medfilt1(data',20)';
                otherwise
                    for j = 1:numCh
                        data(j,:) = lowPassFilterTest(data(j,:),freq/dfreq);
                    end   
             end

            for j = 1:numCh

                %select feature for feature vectors

                switch scoMeth{nthMeth} %subject number = nthDir

                    case 'rmsm'
                    % calculate RMS
                    rmsP = rms(data(j,:), 2); 
                    rmsSco(j) = rmsP;

                    case 'rmsmT3'
                    % calculate RMS average for Top3 / RMS average for all
                    if j == 5 || j == 6
                        rmsSco(j) = rmsT3P;
                    else
                        rmsSco(j) = rmsT3P;
                    end
                    rmsT3P = rms(data(j,:), 2);
                    rmsT3Sco(j) = rmsT3P;     

                    case 'var'
                       varS(j)=var(data(j,:));                

                    case 'mll'
                    % calculate lineLength
                    ll(j) = lineLength(data(j,:));

                    case {'tmmSco', 'smmSco','mspSco'}

                        %peakFinding
                        [pks,locs] = findpeaks(data(j,:),'minpeakdistance',20,'MINPEAKHEIGHT',5);
                        if numel(locs) ~= 0
                            mLocs = max(locs);
                        else
                            mLocs = 0.1;
                        end
                        [npks,nlocs] = findpeaks(mLocs-data(j,:),'minpeakdistance',10,'MINPEAKHEIGHT',5);
                        npks = mLocs-npks;
                        % %plot
                        %figure;plot(x,seg(j,:)), hold on; plot(x(locs),pks+0.05,'k^','markerfacecolor',[0 1 1]); hold on; plot(x(nlocs),npks+0.05,'k^','markerfacecolor',[1 0 0]), hold off

                        %check if peaks founds
                        if numel(nlocs) == 0 
                            switch scoMeth{nthMeth}
                                case 'smmSco'
                                    % no peaks so zero score
                                    mmpSco(j) = 0;
                                case 'tmmSco'
                                    % no peaks so zero score
                                    tmpSco(j) = 0;
                                case 'mspSco'
                                    mspSco(j) = 0;
                            end
                        else

                            switch scoMeth{nthMeth}
                                case 'smmSco'
                                    % calculate sum max/min score
                                    mmpSco(j) = size(pks,2);
                                case 'tmmSco'
                                    % calculate thresholded max/min score
                                    tPks = pks(pks>100);
                                    tmpSco(j) = size(tPks,2);
                                case 'mspSco'
                            
                                    % calculate slopeScore, sS, parameters
                                    a = [pks',locs'];
                                    b = [npks',nlocs'];
                                    sS = vertcat(a,b);
                                    sS = sortrows(sS,2);

                                    % check if slope calcuable 
                                    if numel(sS) > 1
                                        % calculate slopeScore, sS
                                        sS2 = circshift(sS,-1); 
                                        sSS = abs((sS2(:,1) - sS(:,1)) ./ (sS2(:,2) - sS(:,2)));% calc dY / dX
                                        mspSco(j) = mean(sSS);% calculate mean slope score
                                        %sspSco(j) = sum(sSS);% calculate sum slope score
                                    else
                                        mspSco(j) = 0;
                                        %sspSco(j) = 0;
                                    end
                            end
                        end

                    case 'mAs'
                        % calculate mean amplitude score
                        n = 100; %number of max values to find
                        [aa,ix] = sort(data(j,:),'descend');
                        aa = aa(1:n);
                        % ix = ix(1:n);
                        aspSco(j) = sum(aa); 
                end

            end

            % Mean score at each 1sec clip
            switch scoMeth{nthMeth} %subject number = nthDir
                    case 'rmsm'
                        IIlfS(i) = mean(rmsSco); 
                    case 'rmsmT3'
                        %find top 3 values
                        A = rmsT3Sco;
                        Asorted = sort(A(:),'descend');
                        IIlfS(i)=mean(Asorted(1:7))/mean(Asorted(8:end));
                    case 'rmsIZ'           
                        rms5 = rms(data(5,:), 2); 
                        rms6 = rms(data(6,:), 2);
                        rms9 = rms(data(9,:), 2);
                        rms10 = rms(data(10,:), 2);

                        IIlfS(i)=rms5+rms6+rms9+rms10;   
                    case 'phS'
                        IIlfS(i)=sum(data(5,:)-data(6,:));
                       % IIlfS(i)=sum(abs(data(5,:))+abs(data(6,:)));

                    case 'var'
    %                     A = varS;
    %                     Asorted = sort(A(:),'descend');
    %                     IIlfS(i)=mean(Asorted(1:5))/mean(Asorted(5:end));                       
                          IIlfS(i)=varS(5)+varS(6)+varS(9)+varS(10);                   

                    case 'mll'
                        IIlfS(i) = mean(ll);
                    case 'smmSco'
                        A = mmpSco;
                        Asorted = sort(A(:),'descend');
                        IIlfS(i)=mean(Asorted(1:7))/mean(Asorted(8:end));  

    %                    IIlfS(i) = mean(mmpSco); % max/min score
                    case 'tmmSco'
                        A = tmpSco;
                        Asorted = sort(A(:),'descend');
                        IIlfS(i)=mean(Asorted(1:7))/mean(Asorted(8:end));
                    case 'mspSco'
                        IIlfS(i) = mean(mspSco);
            end



        end

        % combine ictal and interictal
        IIIlfS = [IlfS,IIlfS];

        %normalize and save    
        clipScores{nthDir}(1,:) = IIIlfS; % low freq score in row 1 of clipScores
        clipScores{nthDir}(2,:) = IIIlfS./max(IIIlfS); % normalized low freq score in row 1 of clipScores
        clipScores{nthDir}(3,:) = zeros(1,size(clipScores{nthDir},2)); %create placeholder for searchScore
        clipScores{nthDir}(4,:) = zeros(1,size(clipScores{nthDir},2)); %create placeholder for searchScore

        figure;plot(IIIlfS);title(strcat(clipDirNames{nthDir},scoMeth(nthMeth)));axis square; set(gcf,'color','w');drawnow;

        lfS(:,nthMeth) = IIIlfS;
        
    end
end

save('P3scores.mat','lfS')


%% SVM
    sS = (1-testSS(7,1:size(lfS,1),1))';
    % Scores vector == lfS
    X = [lfS(:,1),lfS(:,2),sS];
    
    % Create true labels vector
    labels = [ones(nIClips,1);zeros((nIIClips),1)];    
    
    % train svm model
    svmM = fitcsvm(X,labels,'ClassNames',[0 1],'Standardize',true,'KernelFunction','rbf');
    cSvmM = compact(svmM);
    cSvmM = fitSVMPosterior(cSvmM,X,labels); %returns probabilities

    % predict w/SVM
    [predLabel,postProbs] = predict(cSvmM,X);
    qS = postProbs(:,2);    

% Plot results
     figure('Position', [100, 500, 1100, 400]);
     subplot(1,3,1);plot(qS);title(strcat('SVM Post-Prob. Subj_',clipDirNames(nthDir)));axis square;axis tight;set(gcf,'color','w');drawnow
    subplot(1,3,2);
    for i=1:size(predLabel,1)
        if predLabel(i) ~= labels(i)
            colr = 'r';
        else
            if labels(i) == 1
                colr = 'g';
            else
                colr = 'b';
            end
        end
        scatter3(sS(i),lfS(i,1),lfS(i,2),colr);xlabel('Time');ylabel('LF Score'),zlabel('searchScore');title(strcat('SVM Results. Subj_',clipDirNames(nthDir)));hold on
    end
    axis square;set(gcf,'color','w');  drawnow;  


% Plot training ROC curve
% Create true labels vector
labels = [ones(nIClips,1);zeros((nIIClips),1)];    

% Compute ROC for SVM
posclass = 1;% Define positive class
[X,Y,T,AUC] = perfcurve(labels,qS,posclass);% Generate AUC

% Plot ROC Curve
subplot(1,3,3)

plot(X,Y)
xlabel('False positive rate'); ylabel('True positive rate')
title(strcat('ROC, Subj_',num2str(nthDir)))
annotation('textbox',...
    [0.9 0.1 0.1 0.1],...
    'String',{'AUC:',num2str(AUC)},...
    'FontSize',12,...
    'FontName','Arial',...
    'LineStyle','--',...
    'EdgeColor',[0 0 0],...
    'LineWidth',1,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);
axis square; set(gcf,'color','w');drawnow;    
drawnow;


% Generate threshold based off of AUC
[~,~,~,~,OPTROCPT] = perfcurve(labels,IIIlfS,posclass);
optTPR = OPTROCPT(1,2); %optimal TPR
[TPR,~,T,~] = perfcurve(labels,IIIlfS,posclass,'Xcrit','TPR');
ind = find(TPR==optTPR);
searchScore = mean(T(ind));
disp(strcat('Score for Patient 3:',num2str(searchScore)));

% Save Model (overwrites previous model!)
%save(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(nthDir)),'cSvmM');

warning(s)  % restore the warning state


