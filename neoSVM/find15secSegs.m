% Concatenates Ictal seizure segments


%your matlab directory
        matdir='/Volumes/bobo/';
        clipsdir='kaggleClips/';

% %upper limit for concatenated matrix length, only use if memory problems
% cMaxSize=1e7;

%patients={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};

patients={'Patient_6'};


%not useful, just initializing
offset=zeros(1,4);

for ii=1:length(patients) 

    %concatenate ictal
    %find info
    nIctalClips=length(dir(strcat(matdir,clipsdir,patients{ii},'/*_ictal_s*.mat')));
    [data,freq,latency,channels]=sParLoad(strcat(matdir,clipsdir,patients{ii},filesep,patients{ii},'_ictal_segment_1'));
    nSamplesPerSegment=size(data,2);
    nChannels=size(data,1);
    nthSeizure=0;
    D=zeros(nChannels,nIctalClips*nSamplesPerSegment);

    latSegs = zeros(nIctalClips,1);
    for i=1:nIctalClips
        [data,freq,latency]=sParLoad(strcat(matdir,clipsdir,patients{ii},filesep,patients{ii},'_ictal_segment_',num2str(i)));
        
        disp(latency)
        if latency <=15
            latSegs(i,1) = 1;
        end
    end

    figure; plot(latSegs); 
end