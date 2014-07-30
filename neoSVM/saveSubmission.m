%% save
submissionTable = vertcat(predsSVM{:});
writetable(submissionTable);


disp('done')