%% Make submission from previously generated testSegment searchScores
% add p3 SVM

% your matlab & clips directory
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
nClips=length(clipDirNames);

%check existence of resultsCell
size(resultsCell)
%predsSVM=cell(1,nClips);
%scoresSVM=cell(1,nClips);

tic

s = warning('off','all'); % turn all warnings off
if true,
    disp('Warnings off while running this code!') %lots of warnings from peakPoints not finding a peak
    warning('Here is a warning that doesn''t have an id.')
end

%% use SVM for patient 3
for ii=7:7 
    
    % load appropriate SVM classifier
    load(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(ii)),'cSvmM');    
    
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    testClips=dir(strcat(clipLoc,filesep,'*_test_segment*.mat'));
    nTests=length(testClips);
    scoMeth = {'mspSco','phS'};

    %data stores
    TlfS = zeros(2,nTests); %ictal low freq score    
    seizure=zeros(nTests,1);
    clip = cell(nTests,1);
    
    for nthMeth=1:length(scoMeth);
        
        %initialize prescores
        mspSco = zeros(numCh,1); 
        phS = zeros(numCh,1); %phS = phase score
%             sspSco = zeros(numCh,1);
%             aspSco = zeros(numCh,1);
%             mmpSco = zeros(numCh,1);
%             tmpSco = zeros(numCh,1);
%            rmsSco = zeros(numCh,1); %RMS pre-score    
%            rmsT3Sco = zeros(numCh,1); %RMS pre-score    
%            varS = zeros(numCh,1);
%            rmsIZ = zeros(numCh,1); %phS = ictalZone RMS

        disp(['Finding ' scoMeth{nthMeth} ' searchScore for: ' clipLoc])

        % low Freq scores for ictal clips
        for i=1:nTests
            %create clip name
            clip{i,1}=strcat(clipDirNames{ii},'_test_segment_',num2str(i),'.mat');

            [data,~,freq,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_test_segment_',num2str(i),'.mat'));
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
                        TlfS(nthMeth,i) = mean(rmsSco); 
                    case 'rmsmT3'
                        %find top 7 values
                        A = rmsT3Sco;
                        Asorted = sort(A(:),'descend');
                        TlfS(nthMeth,i)=mean(Asorted(1:7))/mean(Asorted(8:end));
                    case 'rmsIZ'           
                        rms5 = rms(data(5,:), 2); 
                        rms6 = rms(data(6,:), 2); 
                        rms9 = rms(data(6,:), 2); 
                        rms10 = rms(data(6,:), 2); 

                        TlfS(nthMeth,i)=rms5+rms6+rms9+rms10;                     
                    case 'phS'
                        TlfS(nthMeth,i)=sum((data(5,:))-data(6,:));
                         %IlfS(i)=sum(abs(data(5,:))+abs(data(6,:)));

                    case 'mll'
                        TlfS(nthMeth,i) = mean(ll);
                    case 'var'
                        %find top 5 values
    %                     A = varS;
    %                     Asorted = sort(A(:),'descend');
    %                     IlfS(i)=mean(Asorted(1:5))/mean(Asorted(5:end));                   
                          TlfS(nthMeth,i)=varS(5)+varS(6)+varS(9)+varS(10);                   

                    case 'smmSco'

                        A = mmpSco;
                        Asorted = sort(A(:),'descend');
                        TlfS(nthMeth,i)=mean(Asorted(1:7))/mean(Asorted(8:end));                    
                        %IlfS(i) = mean(mmpSco); % max/min score
                    case 'tmmSco'           
                        A = tmpSco;
                        Asorted = sort(A(:),'descend');
                        TlfS(nthMeth,i)=mean(Asorted(1:7))/mean(Asorted(8:end)); 

                    case 'mspSco'
                        TlfS(nthMeth,i) = mean(mspSco);

            end
        end
    end
end   

%% join scores

sS7 = 1-resultsCell{1, 7}.seizure(:,7);
X = [TlfS(1,:);TlfS(2,:);sS7']';

% predict w/SVM
[predLabel,postProbs] = predict(cSvmM,X);
seizure = postProbs(:,2);   
    
%assume early probability is same as seizure
early=seizure;
predsSVM{ii} = table(clip,seizure,early);
scoresSVM{ii} = X;

save('sSVMResults20140726_try2.mat','predsSVM','scoresSVM');
toc 
warning(s)  % restore the warning state
disp('done')




