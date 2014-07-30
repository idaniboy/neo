% define ceiling from cScores

ceilCS = zeros(1,size(clipScores,1));
for i=1:size(clipScores,1)
    meanCS = mean(clipScores{i}(:,1));
    stdCS = std(clipScores{i}(:,1));
    ceilCS(i) = meanCS + 2*stdCS;
end
    
ceilCS