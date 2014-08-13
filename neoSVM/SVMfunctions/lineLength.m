function dist = lineLength(y)
%x and y are vectors containing the coordinates of the points 
%     (in order) between the start and end points
x=1:size(y,2);
dx = diff(x);  %incremental difference between x coordinates
dy = diff(y);  %incremental difference between y coordinates
dL = sqrt(dx.^2+dy.^2);  %length of each segment
dist = sum(dL);    %total line length
end