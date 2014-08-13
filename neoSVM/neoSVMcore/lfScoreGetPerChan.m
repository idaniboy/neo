function prelfS = lfScoreGet(data,meth_ch,numCh,freq,clipLoc)

global matdir
global clipsdir
global clipDirNames
global clipLoc
global scoMeth
global nIClips
global nIIClips
global nTSClips

%data stores
prelfS = []; %ictal low freq score    

%initialize prescores
preSco = zeros(numCh,1);     

%parse meth_ch for per Channel operations
[m_c] = strsplit(meth_ch, ' ');
meth = m_c{1};
if numel(m_c)>1
    chans = strsplit(m_c{2},',');
    if numel(chans)~=0;
        numCh=numel(chans);
    end
    chanList = zeros(1,numel(chans));
    for i = 1:numel(chans)
        chanList(1,i) = str2num(chans{i});
    end    
else
    chanList = 1:numCh;    
end

 
for i = 1:numCh
    
    j= chanList(1,i);
    
    switch freq
        case 5000
            %downsample if freq is 5000
            dfreq = 5;
            Ds = zeros(size(data,1),size(data,2)/5); %assumes SR 5000
            Ds = decimate(data(j,:),dfreq);
    end

    %Pre-process
    switch meth
        case 'phS'
                Ds = medfilt1(Ds',20)';
        otherwise
                Ds = lowPassFilterTest(Ds,freq/dfreq);
    end

    
    %select feature for feature vectors
    switch meth%subject number = nthDir

        case 'rmsm'
        % calculate RMS
        rmsP = rms(Ds, 2); 
        preSco(j) = rmsP;

        case 'rmsmT3'
        % calculate RMS average for Top3 / RMS average for all
        rmsT3P = rms(Ds, 2); 
        if j == 5 || j == 6
            preSco(j) = rmsT3P;
        else
            preSco(j) = rmsT3P;
        end
        preSco(j) = rmsT3P;     

        case 'mll'
            % calculate lineLength
            ll(j) = lineLength(Ds);

        case 'var'
           preSco(j)=var(Ds);

        case {'tmmSco','smmSco','mspSco','mspSco0506','mspSco1011'}

            %peakFinding
            [pks,locs] = findpeaks(Ds,'minpeakdistance',80,'MINPEAKHEIGHT',5);
            if numel(locs) ~= 0
                mLocs = max(locs);
            else
                mLocs = 0.1;
            end
            [npks,nlocs] = findpeaks(mLocs-Ds,'minpeakdistance',80,'MINPEAKHEIGHT',5);
            npks = mLocs-npks;
            % %plot
              %     figure;plot(Ds), hold on; plot(locs,pks+0.05,'k^','markerfacecolor',[0 1 1]); hold on; plot(nlocs,npks+0.05,'k^','markerfacecolor',[1 0 0]), hold off
            %check if peaks founds
            if numel(nlocs) == 0 
                switch meth
                    case 'smmSco'
                        % no peaks so zero score
                        preSco(j) = 0;
                    case 'tmmSco'
                        % no peaks so zero score
                        preSco(j) = 0;
                    case {'mspSco0506','mspSco1011','mspSco'}
                        preSco(j) = 0;
                        %sspSco(j) = 0;
                end
            else

                switch meth
                    case 'smmSco'
                        % calculate sum max/min score
                        preSco(j) = size(pks,2);
                    case 'tmmSco'
                        % calculate thresholded max/min score
                        tPks = pks(pks>100);
                        preSco(j) = size(tPks,2);
                    case {'mspSco','mspSco0506','mspSco1011'}

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
                            preSco(j) = mean(sSS);% calculate mean slope score
                            %sspSco(j) = sum(sSS);% calculate sum slope score
                        else
                            preSco(j) = 0;
                            %sspSco(j) = 0;
                        end
                end
            end

        case 'mAs'
            % calculate mean amplitude score
            n = 100; %number of max values to find
            [aa,ix] = sort(Ds,'descend');
            aa = aa(1:n);
            % ix = ix(1:n);
            preSco(j) = sum(aa); 
    end

end

% Mean score at each 1sec clip
switch meth%subject number = nthDir
        case 'rmsm'
            preSco(preSco==0) = []; %remove zero elements
            prelfS = mean(preSco); 
        case 'rmsmT3'
            %find top 7 values
            A = preSco;
            Asorted = sort(A(:),'descend');
            prelfS=mean(Asorted(1:7))/mean(Asorted(8:end));
        case 'rmsIZ'           
            rms5 = rms(data(5,:), 2); 
            rms6 = rms(data(6,:), 2); 
            rms9 = rms(data(6,:), 2); 
            rms10 = rms(data(6,:), 2); 

            prelfS=rms5+rms6+rms9+rms10;                     
        case 'phS'
            prelfS=abs(sum((data(5,:))-data(6,:)));
             %IlfS(i)=sum(abs(data(5,:))+abs(data(6,:)));

        case 'mll'
            preSco(preSco==0) = []; %remove zero elements
            prelfS = mean(ll);
        case 'var'
            %find top 5 values
%                     A = varS;
%                     Asorted = sort(A(:),'descend');
%                     IlfS(i)=mean(Asorted(1:5))/mean(Asorted(5:end));                   
              prelfS=preSco(5)+preSco(6)+preSco(9)+preSco(10);                   

        case 'smmSco'
            A = preSco;
            Asorted = sort(A(:),'descend');
            prelfS=mean(Asorted(1:7))/mean(Asorted(8:end));                    
            %IlfS(i) = mean(mmpSco); % max/min score
        case 'tmmSco'           
            A = preSco;
            Asorted = sort(A(:),'descend');
            prelfS=mean(Asorted(1:7))/mean(Asorted(8:end)); 

        case 'mspSco0506'
            preSco(preSco==0) = []; %remove zero elements
            prelfS = mean(preSco(5:6));

        case 'mspSco1011'
            preSco(preSco==0) = []; %remove zero elements
            prelfS = mean(preSco(10:11));
            
        case 'mspSco'
            preSco(preSco==0) = []; %remove zero elements
            prelfS = mean(preSco);

end
