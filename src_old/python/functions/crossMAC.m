function MAC = crossMAC(modeShape1,modeShape2,MODE,FREQ)
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

% Plot MAC
figure
barMAC = bar3(MAC);
for k = 1:length(barMAC)   
    zdata = barMAC(k).ZData;   
    barMAC(k).CData = zdata;
    barMAC(k).FaceColor = 'interp';
end
colormap(jet);
colorbar
if MODE==1
    title('MAC - Numerical compared to SSI')
elseif MODE==2
    title('MAC - Numerical compared to ERA')
else
    title('MAC - Numerical compared to FDD')
end
xlabel('OMA frequencies [Hz]')
ylabel('Numerical frequencies [Hz]')
xticks([1,2,3,4,5])
xticklabels(string(FREQ(:,1)'))
yticks([1,2,3,4,5])
yticklabels(string(FREQ(:,2)'))
box on
end