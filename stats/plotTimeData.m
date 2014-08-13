
%best to do one at a time for now, fill in info below

%your matlab directory
matdir='/Volumes/bobo/';
%matdir='C:\Users\paul\Documents\MATLAB\';
%clips dir(relative to matlab, leave off^ that slash, as below. these get mushed together)
clipsdir='kaggleClips/';
%requires data to be in separate folders per patient with _ in folder name
a=dir(strcat(matdir,clipsdir,'*_*'));
clipDirNames={'Patient_3'};%{a.name};%{'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;

for nthDir=1:length(clipDirNames);
    clear save_hashes
    clipLoc=strcat(matdir,clipsdir,clipDirNames{nthDir});
    %load concatenated_ filenames
    nSeizures=length(dir(strcat(clipLoc,filesep,'*_ictal_concatenated_*.mat')));
%         %used to track which seizure within a folder is currently active
%         aSHIndex=1;
    H=[];
    for nthSeizure=1:nSeizures
        figure 

        [data,freq,~,channels]=sParLoad(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_concatenated_',num2str(nthSeizure),'.mat'));
        freq = ceil(freq);

        scaleSpace = 1000;
        colors=cellstr(['black']);
        count=0;

        %colors for each electrode. The order and length must match the electrode
        %list.
        for i=1:size(fieldnames(channels),1)
            colors{i}='black';
        end

        %find data for each channel
        for c = 1:size(fieldnames(channels),1)
            if c == 5 || c == 6
                colorEEG = 'red';
            else
                colorEEG = 'black';
            end
            %plots the data for this electrode at the appropriate y
            %axis location.
            plotLines = plot(data(c,:) + c*scaleSpace, colorEEG);
            
            hold on

        end     
        %% This section fixes the axes labels and ranges 

        % Y axis - uses electrode list as labels 
        maxY = scaleSpace*(length(fieldnames(channels)) + 1);
        ylim([0, maxY])
        set(gca,'YTickLabel',fieldnames(channels))
        set(gca, 'YTick', scaleSpace:scaleSpace:(maxY - scaleSpace) )

    end

    %% Plot interictal data
   nIIClips=length(dir(strcat(clipLoc,filesep,'interictal_concatenated/','*_interictal_concatenated_*.mat')));
%    nIIClips=length(dir('*_interictal_concatenated_*.mat'));

    for nthSeizure=1:nIIClips
        figure 

        [data,freq,~,channels]=sParLoad(strcat(clipLoc,filesep,'interictal_concatenated/',clipDirNames{nthDir},'_interictal_concatenated_',num2str(nthSeizure),'.mat'));
        freq = ceil(freq);

        scaleSpace = 1000;
        colors=cellstr(['black']);
        count=0;

        %colors for each electrode. The order and length must match the electrode
        %list.
        for i=1:size(fieldnames(channels),1)
            colors{i}='black';
        end

        %find data for each channel
        for c = 1:size(fieldnames(channels),1)

            %plots the data for this electrode at the appropriate y
            %axis location.
            if c==5||c==6
                colorEEG = 'red';
            else
                colorEEG = 'black';
            end
            plotLines = plot(data(c,:) + c*scaleSpace, colorEEG);
            
            hold on

        end     
        %% This section fixes the axes labels and ranges 

        % Y axis - uses electrode list as labels 
        maxY = scaleSpace*(length(fieldnames(channels)) + 1);
        ylim([0, maxY])
        set(gca,'YTickLabel',fieldnames(channels))
        set(gca, 'YTick', scaleSpace:scaleSpace:(maxY - scaleSpace) )
        
        drawnow
        pause
        close(gcf)
        
    end
    
end


