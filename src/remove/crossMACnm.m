function MAC = crossMACnm(modeShape1,modeShape2)
% This Matlab script is to calculate and plot cross Modal Assurance Criterion
% (MAC) based on the given mode shapes
% Author: Liangliang Cheng, Postdoctoral researcher, Ghent University
% Email: liangliang.cheng@ugent.be
% Data: 30/11/2020

% Note: madeShape is a matrix with rows (number of DOFs) and columns (number of modes)
% Example: load modeShape1.mat
%          load modeShape2.mat
%          MAC = crossMAC(modeShape1,modeShape2);


[~, numMode1]= size(modeShape1);
[~, numMode2]= size(modeShape2);

MAC = zeros(numMode1,numMode2,'double'); % Initialize MAC matrix
for mode1 = 1:numMode1    
    for mode2 = 1:numMode2
        MAC(mode1,mode2) = abs(modeShape1(:,mode1)'*modeShape2(:,mode2))/ ...
                           sqrt(abs((modeShape1(:,mode1)'*modeShape1(:,mode1))*(modeShape2(:,mode2)'*modeShape2(:,mode2))));                          
    end   
end
