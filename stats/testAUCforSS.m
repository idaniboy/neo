% test searchScore AUC

% Training data
scores=[sS]';

% Create true labels vector
labels = [ones(nIClips,1);zeros((nIIClips),1)];    

%% Find Threshold
% Define positive class
posclass = 1;

% Generate AUC
[~,~,~,~,OPTROCPT] = perfcurve(labels,scores,posclass);
optTPR = OPTROCPT(1,2); %optimal TPR

[TPR,~,T,~] = perfcurve(labels,scores,posclass,'Xcrit','TPR');
ind = find(TPR==optTPR);
searchScoreThresh = mean(T(ind));



%% Compute ROC 
% Define positive class
posclass = 1;

% Generate AUC
[X,Y,T,AUC] = perfcurve(labels,scores,posclass);

% Plot ROC Curve
figure;
plot(X,Y)
xlabel('False positive rate'); ylabel('True positive rate')
title('ROC, using mean landmarks per top 200 matches, basis X-axis')
annotation('textbox',...
    [0.3 0.65 0.2 0.1],...
    'String',{'AUC:',num2str(AUC)},...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','--',...
    'EdgeColor',[0 0 0],...
    'LineWidth',1,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);
axis square; set(gcf,'color','w');

% Generate gamma score
[~,~,~,~,OPTROCPT] = perfcurve(labels,scores,posclass);
[TPR,~,T,~] = perfcurve(labels,scores,posclass,'Xcrit','TPR');
gammaScore(nthDir) = mean(T(TPR==OPTROCPT(1,2)));
disp(strcat('Gamma score: ',num2str(gammaScore(nthDir))));
toc





