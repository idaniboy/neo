% 1st section trains SVM using searchScore and lowFreqScore
% requires clipScores
% test run using pre-calculated searchScore values from prior run
% for normalized searchScores and NON-normalized lowFreqScores

%size(clipScores)

% your matlab & clips directory
matdir='/Volumes/bobo/';
clipsdir='kaggleClips/';
clipDirNames={'Dog_1','Dog_2','Dog_3','Dog_4','Patient_1','Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8'};%,'Patient_6}%;
nClips=length(clipDirNames);

% define ceiling of clipScores for normalization
% ceilCS = zeros(1,size(clipScores,1));
% for i=1:size(clipScores,1)
%     meanCS = mean(clipScores{i}(1,:));
%     stdCS = std(clipScores{i}(1,:));
%     maxCS = max(clipScores{i}(1,:));
%     %ceilCS(i) = meanCS + 2*stdCS;
%     ceilCS(i) = maxCS;
%     cS = clipScores{i}(1,:);
%     %cS(cS>ceilCS(i)) = ceilCS(i);
%     cS = cS/ceilCS(i);
%     clipScores{i}(1,:) = cS;
% end

% normalize clipScores

%% Train SVM
% for real data, still uses same parameters
for ii=7:7

    % init vars for test subject
    clipLoc=strcat(matdir,clipsdir,clipDirNames{ii});
    nIClips=length(dir(strcat(clipLoc,filesep,'*_ictal_segment_*.mat')));
    nIIClips=length(dir(strcat(clipLoc,filesep,'*_interictal_segment_*.mat')));
    nIIIClips = nIClips + nIIClips;
    
    % load feature vectors
%    lfS = clipScoresP3{ii}(1,:);
    %sS = clipScoresSS{ii}(2,:);
    %lfS_phS = p3Scores(1,:);
%     %lfS_smpSco = p3Scores(2,:);
%     mpD = clipBpScores{ii}(1,:);
%     mpB = clipBpScores{ii}(1,:);
%     mpG = clipBpScores{ii}(1,:);
    
    % Put feature vectors into X
    switch ii
        case 7
            % Scores vector == lfS
                sS = (1-testSS(7,1:size(lfS,1),1))';
    % Scores vector == lfS
    X = [lfS(:,1),lfS(:,2),sS];
            
        case 12
            if ii == 12
                %sS = 1-clipScoresSS{ii}(2,:);
                X=[mpB;mpG;lfS]';
            else
                X=[mpD;lfS;sS]';
            end
    end
    
    % Create true labels vector
    labels = [ones(nIClips,1);zeros((nIIClips),1)];    
       
    % train svm model
    svmM = fitcsvm(X,labels,'ClassNames',[0 1],'Standardize',true,'KernelFunction','rbf');
    cSvmM = compact(svmM);
    cSvmM = fitSVMPosterior(cSvmM,X,labels); %returns probabilities

    % predict w/SVM
    [predLabel,postProbs] = predict(cSvmM,X);
    qS = postProbs(:,2);    

%     % Plot results
     figure('Position', [100, 500, 1100, 400]);
     subplot(1,3,1);plot(qS);title(strcat('SVM Post-Prob. Subj_',num2str(ii)));axis square;axis tight;set(gcf,'color','w');drawnow
%     subplot(1,3,2);
%     for i=1:size(predLabel,1)
%         if predLabel(i) ~= labels(i)
%             colr = 'r';
%         else
%             if labels(i) == 1
%                 colr = 'g';
%             else
%                 colr = 'b';
%             end
%         end
%         scatter3(i,lfS(i),sS(i),colr);xlabel('Time');ylabel('LF Score'),zlabel('searchScore');title(strcat('SVM Results. Subj_',num2str(ii)));hold on
%     end
%     axis square;set(gcf,'color','w');  drawnow;  
    
    % Compute ROC for SVM
    posclass = 1;% Define positive class
    [X,Y,T,AUC] = perfcurve(labels,qS,posclass);% Generate AUC

    % Plot ROC Curve
    subplot(1,3,3);
    plot(X,Y);    drawnow;
    xlabel('False positive rate'); ylabel('True positive rate')
    title(strcat('ROC, Subj_',num2str(ii)))
    drawnow;
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
    save(strcat('/Users/idaniboy/Documents/MATLAB/svmTrained/svmM_Subject_',num2str(ii)),'cSvmM');
    drawnow;
end
