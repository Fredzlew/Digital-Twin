% parameters
clc; clear; close all;
addpath(genpath('..\..\data'))
% Loading modal parameters from OMA 
promptt = "High damping or no damping? (1 = High and 2 = no damp): ";
xx = input(promptt);
if xx == 1
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par\SSIfreq_5_2_1.npy');
    SSIomega = SSIFreq * 2 * pi;
    SSIphi = readNPY('..\..\data\experimental_data\Modal_par\SSImodes_5_2_1.npy');
elseif xx == 2
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par\SSIfreq_no_damp.npy');
    SSIomega = SSIFreq * 2 * pi;
    SSIphi = readNPY('..\..\data\experimental_data\Modal_par\SSImodes_no_damp.npy');
end



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

% damping different





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
            legend('SSI','Location','northwest')
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
OMAfreq = SSIFreq;


%% Download the file to plot in latex
%{
if xx == 1
    T = array2table([num2cell(x'),num2cell(phi),num2cell([0;OMAfreq])]);
    T.Properties.VariableNames(1:7) = {'height','OMAphi1','OMAphi2','OMAphi3','OMAphi4','OMAphi5','OMAfreq'};
    writetable(T,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap5_SSIcov_on_experimental_data_highdamp.csv','delimiter',';')
elseif xx == 2
    T = array2table([num2cell(x'),num2cell(phi),num2cell([0;OMAfreq])]);
    T.Properties.VariableNames(1:7) = {'height','OMAphi1','OMAphi2','OMAphi3','OMAphi4','OMAphi5','OMAfreq'};
    writetable(T,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap5_SSIcov_on_experimental_data_nodamp.csv','delimiter',';')
end
%}