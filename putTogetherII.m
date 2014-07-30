% Concatenates interictal seizure segments


%your matlab directory
        matdir='/Volumes/bobo/';
        clipsdir='kaggleClips/';

% %upper limit for concatenated matrix length, only use if memory problems
% cMaxSize=1e7;

%patients={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};
patients={'Patient_3'};
%patients={'Patient_8'};

for ii=1:length(patients) 

    %concatenate ictal
    %find info
    nIIctalClips=length(dir(strcat(matdir,clipsdir,patients{ii},'/*_interictal_s*.mat')));
    [data,freq,latency,channels]=sParLoad(strcat(matdir,clipsdir,patients{ii},filesep,patients{ii},'_interictal_segment_1.mat'));
    nSamplesPerSegment=size(data,2);
    nChannels=size(data,1);
    nthSeizure=0;
    
    sampPerSeg = 10;
    j=0;
    D=zeros(nChannels,sampPerSeg*nSamplesPerSegment);

    for i=1:nIIctalClips
        [data,freq,latency]=sParLoad(strcat(matdir,clipsdir,patients{ii},filesep,patients{ii},'_interictal_segment_',num2str(i),'.mat'));
        
        %write current seizure data into D starting at D(:,1)
        D(:,((i-j*sampPerSeg)-1)*nSamplesPerSegment+1:(i-j*sampPerSeg)*nSamplesPerSegment)=data;
  
        if i ==sampPerSeg*(j+1)
            save(strcat(patients{ii},'_interictal_concatenated_',num2str(j+1),'.mat'),'D','channels','freq');
            j = j + 1;
            D=zeros(nChannels,sampPerSeg*nSamplesPerSegment);    
        end
    end
    
    
end