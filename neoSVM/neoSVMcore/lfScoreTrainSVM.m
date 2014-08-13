function X = lfScoreTrainSVM(subj,subjEnd)
%% creates feature vectors from low frequency features of seizure data
%% then trains SVM model
% derived from lfSDforClips
% chooses feature based off performance on training data
% plots AUC curves for each subject
% set subj to 1 and subjEnd to 12 to run through all subjects
% set subj == subjEnd to run just one subject

% Init vars
global matdir
global clipsdir
global clipDirNames
global clipLoc
global scoMeth
global nIClips
global nIIClips
global nTSClips

s = warning('off','all'); % turn all warnings off
if true,
    disp('Warnings off while running this code!') %lots of warnings from peakPoints not finding a peak
    warning('Here is a warning that doesn''t have an id.')
end

for nthDir=subj:subjEnd;  

    %data store for lfS
    lfS = zeros(length(scoMeth{nthDir}),nIClips+nIIClips);
    
    for nthMeth=1:length(scoMeth{nthDir});

        %data stores
        IlfS = zeros(1,nIClips); %ictal low freq score
        IIlfS = zeros(1,nIIClips); %interictal low freq score
        disp(['Finding ' scoMeth{nthDir}{nthMeth} ' searchScore for: ' clipLoc])

        % Scores for ictal clips
        for i = 1:nIClips;
            [data,~,freq,channels]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_ictal_segment_',num2str(i),'.mat'));
            numCh = length(fieldnames(channels));
            IlfS(1,i) = lfScoreGetPerChan(data,scoMeth{nthDir}{nthMeth},numCh,freq,clipLoc);
        end

        
        % Scores for interictal clips
        for i = 1:nIIClips;
            [data,~,freq,channels]=sParLoad_oldV(strcat(clipLoc,filesep,clipDirNames{nthDir},'_interictal_segment_',num2str(i),'.mat'));
            numCh = length(fieldnames(channels));
            IIlfS(1,i) = lfScoreGetPerChan(data,scoMeth{nthDir}{nthMeth},numCh,freq,clipLoc);
        end

        % combine ictal and interictal
        IIIlfS = [IlfS,IIlfS];
        %lfS(:,nthMeth) = IIIlfS;
        lfS(nthMeth,:) = IIIlfS;
        
        figure;plot(IIIlfS);title(strcat(clipDirNames{nthDir},scoMeth{nthDir}(nthMeth)));axis square; set(gcf,'color','w');drawnow;
        
    end
end

save(strcat('P',num2str(subj),'_scores.mat'),'lfS')

%% SVM for seizure or not
%If using searchScore
%sS = (1-testSS(10,1:size(lfS,2),1))';
switch nthDir
    case 6
        % code for making ratios from first 3 feature pairs for P2
        Xtemp = lfS; %X is longest vertically
        X = zeros(3,size(lfS,2));
        X(1,:) = Xtemp(1,:)./Xtemp(2,:);
        X(2,:) = Xtemp(3,:)./Xtemp(4,:);
        X(3,:) = Xtemp(5,:)./Xtemp(6,:);
        X = X';
    otherwise
        X = lfS';        
end

% Create true labels vector
labels = [ones(nIClips,1);zeros((nIIClips),1)];    

% train svm model
%SVMModel = fitcsvm(X,labels,'KernelFunction','rbf','OutlierFraction',0.05,'boxconstraint',0.01);
SVMModel = fitcsvm(X,labels,'KernelFunction','rbf','OutlierFraction',0.01,'boxconstraint',10,'crossval','on','Holdout',0.2,'Cost',[0,1;2,0]);
disp('kfoldLoss');
kfoldLoss(SVMModel);
cSvmM = compact(SVMModel);
cSvmM = fitSVMPosterior(cSvmM,X,labels); %returns probabilities

% Save Model (overwrites previous model!)
save(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(nthDir)),'cSvmM','SVMModel');

% predict w/SVM
[predLabel,postProbs] = predict(cSvmM,X);
qS = postProbs(:,2);    

% Plot results
figure('Position', [100, 500, 1100, 400]);
subplot(1,3,1);plot(qS);title(strcat('SVM Post-Prob. Subj_',clipDirNames(nthDir)));axis square;axis tight;set(gcf,'color','w');drawnow
% subplot(1,3,2);
% for i=1:size(predLabel,1)
%     if predLabel(i) ~= labels(i)
%         colr = 'r';
%     else
%         if labels(i) == 1
%             colr = 'g';
%         else
%             colr = 'b';
%         end
%     end
%     scatter3(sS(i),lfS(i,1),lfS(i,2),colr);xlabel('searchScore');ylabel('slopeScore'),zlabel('phaseScore');title(strcat('SVM Results. Subj_',clipDirNames(nthDir)));hold on
% end
% axis square;set(gcf,'color','w');  drawnow;  


% Plot training ROC curve
% Create true labels vector
labels = [ones(nIClips,1);zeros((nIIClips),1)];    

% Compute ROC for SVM
posclass = 1;% Define positive class
[Xvals,Y,T,AUC] = perfcurve(labels,qS,posclass);% Generate AUC

% Plot ROC Curve
subplot(1,3,3)

plot(Xvals,Y)
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
% 
warning(s)  % restore the warning state

