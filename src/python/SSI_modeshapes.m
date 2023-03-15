
% parameters
clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'),genpath('python'))
% Loading modal parameters from OMA 
SSIFreq = readNPY('omega.npy');
SSIomega = SSIFreq * 2 * pi;
SSIphi = readNPY('phi.npy');


% dimensions in meters
t = 0.015; % floor height [m]
h = 1*10^-3; % short side of column [m]
b = 30*10^-3; % long side (depth) of column [m]
L = 175*10^-3; % column length in middle [m]
Lb = 168*10^-3; % column length at bottom [m]
Lt = 75*10^-3; % column length on top [m]

% story heights [m] (from ground to mid floor)
H(1) = Lb + t/2;
for i = 2:5
    H(i) = H(i-1) + L + t;
end

% plotting the mode shapes
phi = [zeros(1,length(SSIphi)); SSIphi];
x = [0, H];
fig = figure;
fig.Position=[100 100 1600 700];
for i=1:5
    subplot(1,5,i)
    hold on
    plot(phi(:,i),x,'-m')
    plot(phi(2:end,i),x(2:end),'b.','markersize',30)
    title(['f = ' num2str(SSIFreq(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
    xline(0.0,'--')
    xlim([-1.1,1.1])
    ylim([0,x(end)])
        if i==1
            legend('SSI','Approximation','Location','northwest')
        else
        end
end
sgtitle('SSI','FontSize',20) 
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);