function [MSacc,TOTacc]=modeshapeacc(v1,v2)
% Define matrix with accuracies
Macc = zeros(size(v1,1),size(v1,2));
% Looping over columns and rows
for i=1:size(v1,1)
    for j=1:size(v1,2)
        Macc(i,j) = min(abs(v1(i,j)),abs(v2(i,j)))/max(abs(v1(i,j)),abs(v2(i,j)));
    end
end
% Final accuracy is averaged over all floors for each mode shape
MSacc = mean(Macc);
% Total accuracy, based on all mode shapes
TOTacc = mean(MSacc);
end