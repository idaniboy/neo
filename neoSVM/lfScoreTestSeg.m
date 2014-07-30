function lfS = lfScoreTestSeg(nthDir,clipLoc,clipName)
% % creates feature vectors from low frequency features of seizure data
% development for kaggle seizure contest...
% derived from lfScore
% chooses feature based off performance on training data

% feature vector order sequence; in future automate discovery of this
% using lfSDforClipsPreTesting
fvOrder = {'mll','smmSco','rmsm','tmmSco','mll','mll','mspSco','rmsm','mll','rmsm','rmsm','rmsm'};
    
%load data
[data,~,freq,channels]=sParLoad_oldV(strcat(clipLoc,filesep,clipName));
data = medfilt1(data',10)'; %consider adjust medfilt1 based on frequency in future
freq=ceil(freq);
numCh=length(fieldnames(channels));

%downsample 5000Hz --> 1000Hz
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

%initialize prescores
mspSco = zeros(numCh,1); % 
sspSco = zeros(numCh,1);
aspSco = zeros(numCh,1);
mmpSco = zeros(numCh,1);
tmpSco = zeros(numCh,1);
rmsSco = zeros(numCh,1); %RMS pre-score    

for j = 1:length(fieldnames(channels))

    %select feature for feature vectors
    switch fvOrder{nthDir} %subject number = nthDir

        case 'rmsm'
        % calculate RMS
        rmsP = rms(data(j,:), 2); 
        rmsSco(j) = rmsP;

        case 'mll'
        % calculate lineLength
        ll(j) = lineLength(data(j,:));

        case {'tmmSco','smmSco'}

            %peakFinding
            [pks,locs] = findpeaks(data(j,:),'minpeakdistance',20,'MINPEAKHEIGHT',20);
            if numel(locs) ~= 0
                mLocs = max(locs);
            else
                mLocs = 0.1;
            end
            [npks,nlocs] = findpeaks(mLocs-data(j,:),'minpeakdistance',10,'MINPEAKHEIGHT',20);
            npks = mLocs-npks;
            % %plot
            %figure;plot(x,seg(j,:)), hold on; plot(x(locs),pks+0.05,'k^','markerfacecolor',[0 1 1]); hold on; plot(x(nlocs),npks+0.05,'k^','markerfacecolor',[1 0 0]), hold off

            %check if peaks founds
            if numel(nlocs) == 0 
                switch fvOrder{nthDir}
                    case 'smmSco'
                        % no peaks so zero score
                        mmpSco(j) = 0;
                    case 'tmmSco'
                        % no peaks so zero score
                        tmpSco(j) = 0;

                        % mspSco(j) = 0;
%                         sspSco(j) = 0;
                end
            else

                switch fvOrder{nthDir}
                    case 'smmSco'
                        % calculate sum max/min score
                        mmpSco(j) = size(pks,2);
                    case 'tmmSco'
                        % calculate thresholded max/min score
                        tPks = pks(pks>100);
                        tmpSco(j) = size(tPks,2);
                end
%                         % calculate slopeScore, sS, parameters
%                         a = [pks',locs'];
%                         b = [npks',nlocs'];
%                         sS = vertcat(a,b);
%                         sS = sortrows(sS,2);
% 
%                         % check if slope calcuable 
%                         if numel(sS) > 1
%                             % calculate slopeScore, sS
%                             sS2 = circshift(sS,-1); 
%                             sSS = abs((sS2(:,1) - sS(:,1)) ./ (sS2(:,2) - sS(:,2)));% calc dY / dX
%                             mspSco(j) = mean(sSS);% calculate mean slope score
%                             sspSco(j) = sum(sSS);% calculate sum slope score
%                         else
%                             mspSco(j) = 0;
%                             sspSco(j) = 0;
%                         end
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
switch fvOrder{nthDir} %subject number = nthDir
        case 'rmsm'
            lfS = mean(rmsSco); 
        case 'mll'
            lfS = mean(ll);
        case 'smmSco'
            lfS = mean(mmpSco); % max/min score
        case 'tmmSco'
            lfS = mean(tmpSco); % thresholded max/min score
end



