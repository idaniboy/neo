function SVMlatencyTrain(X,ii)
%% SVM for seizure latency
% requires that lfScoreTrainSVM be run first to generate X (scores from...
% feature vectors)
% ii = subj number

% Init vars
global matdir
global clipsdir
global clipDirNames
global clipLoc
global scoMeth
global nIClips
global nIIClips
global nTSClips

latSegs = zeros(nIClips,1);
for i=1:nIClips
    [~,~,latency]=sParLoad(strcat(matdir,clipsdir,clipDirNames{ii},filesep,clipDirNames{ii},'_ictal_segment_',num2str(i)));
    if latency <=15
        latSegs(i,1) = 1;
    end
end
figure; plot(latSegs); title(strcat('Latency for ',num2str(ii)));

% Create scores vector
latScores = X(1:nIClips,:);

% Create true labels vector
labels = latSegs;    

% train svm model
SVMLatModel = fitcsvm(latScores,labels,'ClassNames',[0 1],'KernelFunction','rbf','OutlierFraction',0.05,'boxconstraint',0.01);
%cSvmM = crossval(SVMModel);
cLatSvmM = compact(SVMLatModel);
cLatSvmM = fitSVMPosterior(cLatSvmM,latScores,labels); %returns probabilities

% Save Model (overwrites previous model!)
save(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Latency_Subject_',num2str(ii)),'cLatSvmM');

% predict w/SVM
[predLabel,postProbs] = predict(cLatSvmM,latScores);
qS = postProbs(:,2);    

% Plot results
figure('Position', [100, 500, 1100, 400]);
subplot(1,3,1);plot(qS);title(strcat('SVM Post-Prob. Subj_',clipDirNames(ii)));axis square;axis tight;set(gcf,'color','w');drawnow
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
%     scatter3(sS(i),lfS(i,1),lfS(i,2),colr);xlabel('searchScore');ylabel('slopeScore'),zlabel('phaseScore');title(strcat('SVM Results. Subj_',clipDirNames(ii)));hold on
% end
% axis square;set(gcf,'color','w');  drawnow;  


% Plot training ROC curve

% Compute ROC for SVM
posclass = 1;% Define positive class
[xVals,Y,T,AUC] = perfcurve(labels,qS,posclass);% Generate AUC

% Plot ROC Curve
subplot(1,3,3)

plot(xVals,Y)
xlabel('False positive rate'); ylabel('True positive rate')
title(strcat('ROC, Subj_',num2str(ii)))
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