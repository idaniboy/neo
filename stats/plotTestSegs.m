
%best to do one at a time for now, fill in info below

%your matlab directory
matdir='/Volumes/bobo/';
%matdir='C:\Users\paul\Documents\MATLAB\';
%clips dir(relative to matlab, leave off^ that slash, as below. these get mushed together)
clipsdir='kaggleClips/';
%requires data to be in separate folders per patient with _ in folder name
a=dir(strcat(matdir,clipsdir,'*_*'));
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
subj = 7;

% parameters
scaleSpace = 1000;
clipLoc=strcat(matdir,clipsdir,clipDirNames{subj});
nTotSegs=length(dir(strcat(clipLoc,filesep,'*_test_segment_*.mat')));
[~,freq,~,channels]=sParLoad(strcat(clipLoc,filesep,clipDirNames{subj},'_test_segment_',num2str(1),'.mat'));
freq = ceil(freq);
numCh = length(fieldnames(channels));
dfreq = 5;

% init figure
figure('units','normalized','position',[0 .9 1 .8])

for nSeg=1:nTotSegs
    xPos = rem(nSeg,10);
    if xPos == 0
        xPos = 10;
    end
    subplot(1,10,xPos)

    [data,~,~,~]=sParLoad(strcat(clipLoc,filesep,clipDirNames{subj},'_test_segment_',num2str(nSeg),'.mat'));

    for j=1:numCh
        %switch freq
            %case 5000
            %downsample if freq is 5000
            Ds = decimate(data(j,:),dfreq);
        %end

        if j == 5 || j == 6
            colorEEG = 'red';
        else
            colorEEG = 'black';
        end
        %plots the data for this electrode at the appropriate y
        %axis location.
        plotLines = plot(Ds + j*scaleSpace, colorEEG);
        hold on            
    end

    % This section fixes the axes labels and ranges 
    % Y axis - uses electrode list as labels 
    if rem(nSeg,10)== 1
        maxY = scaleSpace*(length(fieldnames(channels)) + 1);
        ylim([0, maxY])
        set(gca,'YTickLabel',fieldnames(channels))
        set(gca, 'YTick', scaleSpace:scaleSpace:(maxY - scaleSpace) )
    else
        set(gca,'YTickLabel','')
    end

    probSz = sprintf('%.2f',resultsCell{1, 7}.seizure(nSeg));
    title(strcat(num2str(nSeg),':',num2str(probSz)));
    axis tight
    drawnow

    if xPos==10
        pause
        close all
        figure('units','normalized','position',[0 .9 1 .8])
    end
end




