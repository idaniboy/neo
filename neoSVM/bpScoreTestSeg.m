function [bpD,bpB,bpG] = bpScoreTestSeg(nthDir,clipLoc,clipName)
% % creates feature vectors from low frequency features of seizure data
% development for kaggle seizure contest...
% derived from lfScore
% chooses feature based off performance on training data

    
%load params
[datum,~,freq,channels]=sParLoad_oldV(strcat(clipLoc,filesep,clipName));
freq=ceil(freq);
numCh = length(fieldnames(channels));

switch freq
    case 400
        zB = zeros(1,312); % total length = 1024
        topFreq = 200;        
    case 500
        zB = zeros(1,262); % total length = 1024
        topFreq = 250;
    case 5000
        zB = zeros(1,12);% total length = 1024
        topFreq = 400;
        %downsample freq from 5000 --> 1000;
        dfreq = 5;
        Ds = zeros(size(datum,1),size(datum,2)/5); %assumes SR 5000
        for j=1:size(datum,1)
            Ds(j,:) = decimate(datum(j,:),dfreq);
        end
        datum = Ds;
        freq = 1000;
end
data = zeros(numCh,(size(zB,2)*2+freq));
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
bpD = mean(ppD);
bpB = mean(ppB);        
bpG = mean(ppG);
%mpIV(i) = mean(ppV);
    
    
end
