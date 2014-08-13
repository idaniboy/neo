% % creates feature vectors from low frequency features of seizure data
% development for kaggle seizure contest...
% derived from lfSDforClips
% chooses feature based off performance on training data
% plots AUC curves for each subject
% non-normalized LFS on row 1
% normalizes samples on 0 to 1 scale for row 2

    %your matlab directory
%    matdir='/Users/idaniboy/Documents/MATLAB/';
    %matdir='C:\Users\paul\Documents\MATLAB\';
    %clips dir(relative to matlab, leave off^ that slash, as below. these get mushed together)
    %clipsdir='kaggleShazam\clips\';
%    clipsdir='kaggleClips/ictalSegs/';
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
    fvOrder = {'mll','smmSco','rmsm','tmmSco','mll','mll','var','rmsm','mll','rmsm','rmsm','rmsm'};

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
        
        %data store for lfS
        lfS = zeros(nIClips+nIIClips,1);
                
        %initialize prescores
        mspSco = zeros(numCh,1); % 

        %data stores
        IlfS = zeros(1,nIClips); %ictal low freq score
        IIlfS = zeros(1,nIIClips); %interictal low freq score

        disp(['Finding searchScore for: ' clipLoc])

        % low Freq scores for ictal clips
        for i = 1:nIClips;
            [data,~,~,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_segment_',num2str(i),'.mat'));
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

            data = medfilt1(data',20)';
            lenD = size(data,2);


            for j = 1:numCh
                                %initialize l and m
                l=1:lenD+2;
                m = ones(1,lenD+2);

                clipImg=[];
                clipImg2=[];
                clipImg3=[];
                clipImg4=[];
                
                % Hough transform
                m(2:lenD+1)=ceil(10*data(j,:));
                maxM = max(m);
                minM = abs(min(m));
                m = m + minM + 1;
                S = sparse(l,m,1,lenD+2,minM+maxM+1);
                clipPts = full(S)';
%                
                clipImg = roipoly(clipPts, l, m);
%                 
%                 clipImg = bwmorph(BW,'remove');
%                 clipImg = bwmorph(clipImg,'bridge',100);
%                 figure, imshow(clipImg)
                
                clipImg2 = [zeros(minM+maxM+1,1),clipImg,zeros(minM+maxM+1,1)];
                clipImg3 = clipImg2 + [clipImg,zeros(minM+maxM+1,1),zeros(minM+maxM+1,1)];
                clipImg4 = clipImg3 + [zeros(minM+maxM+1,1),zeros(minM+maxM+1,1),clipImg];
                clipImg4(clipImg4==3)=0;
                clipImg4(clipImg4==2)=1;
%                figure;imagesc(clipImg4);
%                  imdata=bwmorph(clipImg,'skel', inf);
%                  figure;imshow(imdata);
                 
%                 
%                 BW = edge(clipImg4,'canny');
%                 figure;imagesc(BW)

                [H,T,R] = hough(clipImg4);
%                 imshow(H,[],'XData',T,'YData',R,...
%                             'InitialMagnification','fit');hold on
%                 xlabel('\theta'), ylabel('\rho');
%                 axis on, axis normal, hold on;
                P  = houghpeaks(H,100,'threshold',ceil(0.2*max(H(:))));
%                 x = T(P(:,2)); y = R(P(:,1));
%                 figure;plot(x,y,'s','color','white');
                % Find lines and plot them
                lines = houghlines(clipImg4,T,R,P,'FillGap',100,'MinLength',100);
                figure, imshow(clipImg4), hold on
                max_len = 0;
                for k = 1:length(lines)
                   xy = [lines(k).point1; lines(k).point2];
                   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');hold on

                   % Plot beginnings and ends of lines
                   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');hold on
                   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');hold on

                   % Determine the endpoints of the longest line segment
                   len = norm(lines(k).point1 - lines(k).point2);
                   if ( len > max_len)
                      max_len = len;
                      xy_long = xy;
                   end
                end

                % highlight the longest line segment
                plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');hold on
                
                % initialize               
                mxy1 = [lines.point1];
                mxy2 = [lines.point2];   
                mxy1 = reshape(mxy1,[2,size(mxy1,2)/2])';
                mxy2 = reshape(mxy2,[2,size(mxy2,2)/2])';                
                
                %sort individual points of each line in order by time
                for k = 1:size(mxy1,1)
                    swapStore1 = mxy1(k,:);
                    swapStore2 = mxy2(k,:);
                    if mxy1(k,1)>mxy2(k,1)
                        mxy1(k,:) = swapStore2;
                        mxy2(k,:) = swapStore1;
                    end
                end                
                
                %sort rows in order by time
                [sxy1, sI] = sortrows(mxy1,1);
                sxy2 = mxy2(sI,:);

                % find slopes of lines
                xySlp=zeros(1,size(sxy1,1));
                for k = 1:size(sxy1,1)
                   % first determine slopes of lines
                   xySlp(k) = (sxy2(k,2)-sxy1(k,2))/(sxy2(k,1)-sxy1(k,1));
                end
                
                fxy1 = [];
                fxy2 = [];
                k = 1;
                % merge consecutive lines with same slope
                while k <= size(sxy1,1)-1
                    if 1==sign(xySlp(k)*xySlp(k+1))

                        
                        fxy1(k,:) = sxy1(k,:);
                        fxy2(k,:) = sxy2(k+1,:);
                        k = k + 1;
                    else
                        fxy1(k,:) = sxy1(k,:);
                        fxy2(k,:) = sxy2(k,:);                        
                    end
                        k = k + 1;
                end
                
                     sxy1(all(sxy1==0,2),:)=[];
                     sxy2(all(sxy2==0,2),:)=[];
                figure, imshow(clipImg4), hold on
                
                for k = 1:length(sxy1)
                   xy = [fxy1(k,:); fxy2(k,:)];
                   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','magenta');hold on

                   % Plot beginnings and ends of lines
                   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');hold on
                   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');hold on

                   % Determine the endpoints of the longest line segment
                   len = norm(fxy1(k) - fxy2(k));
                   if ( len > max_len)
                      max_len = len;
                      xy_long = xy;
                   end
                   
                   % highlight the longest line segment
                   %plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');
                   pause
                end               
                
            end
           
        end

        % mpIIMG = zeros(1,nIIClips);
        % mpIIHG = zeros(1,nIIClips);

        %initialize prescores
        mspSco = zeros(numCh,1); % 
        
        % bandpower for interictal clips
        for i = 1:nIIClips;
            [data,~,~,~]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_interictal_segment_',num2str(i),'.mat'));
            data = medfilt1(data',10)';

            for j = 1:numCh

                %select feature for feature vectors

               

            end

            % Mean score at each 1sec clip
       


        end

        % combine ictal and interictal
        IIIlfS = [IlfS,IIlfS];

        %normalize and save    
        clipScores{nthDir}(1,:) = IIIlfS; % low freq score in row 1 of clipScores

        figure;plot(IIIlfS);title(strcat(clipDirNames{nthDir}));axis square; set(gcf,'color','w');drawnow;

        lfS(:,1) = IIIlfS;
        
    
end

save('cScores.mat','clipScores')


%% Plot training ROC curve
% Create true labels vector
labels = [ones(nIClips,1);zeros((nIIClips),1)];    

% Compute ROC for SVM
posclass = 1;% Define positive class
[X,Y,T,AUC] = perfcurve(labels,lfS,posclass);% Generate AUC

% Plot ROC Curve
%subplot(1,3,3)
figure
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

warning(s)  % restore the warning state


