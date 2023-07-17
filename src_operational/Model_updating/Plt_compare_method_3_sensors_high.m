%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting the model update %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath(genpath('..\data'),genpath('.\Costfunctions\data_updated_par_3_sensors'),genpath('.\Costfunctions\functions'),genpath('.\Sensitivity_method\data_updated_par_sens_3_sensors'))

% Loading stiffness for all OMA costfunction and numerical model
filename = load('..\data\modelprop.mat'); % omegas from numericla model
Globalstiff = filename.k;

% Costfunctions
dataSSIFreq = load('.\Costfunctions\data_updated_par_3_sensors\SSIfreq_3_sensors_high.mat'); % SSI FREQ
K_SSI_freq = dataSSIFreq.k;

dataSSImode = load('.\Costfunctions\data_updated_par_3_sensors\SSImode_3_sensors_high.mat'); % SSI modes
K_SSI_mode = dataSSImode.k;

dataSSImodemac = load('.\Costfunctions\data_updated_par_3_sensors\SSImode_mac_3_sensors_high.mat'); % SSI modes
K_SSI_mode_mac = dataSSImodemac.k;


dataSSIFreqmode = load('.\Costfunctions\data_updated_par_3_sensors\SSIfreqmode_3_sensors_high.mat'); % SSI FREQ and modes
K_SSI_freq_mode = dataSSIFreqmode.k;

dataSSIfreqmodeEIL = load('.\Costfunctions\data_updated_par_3_sensors\SSIfreqmodeEIL_3_sensors_high.mat'); % SSI modes
K_SSI_freqmodeEIL = dataSSIfreqmodeEIL.k;


% Sensitivity method
dataSSIFreq_sens = load('.\Sensitivity_method\data_updated_par_sens_3_sensors\Eigenvalue_residual_3_sensors_high.mat'); % SSI FREQ
K_SSI_freq_sens = dataSSIFreq_sens.knew';

dataSSImode_sens = load('.\Sensitivity_method\data_updated_par_sens_3_sensors\Mode_shape_residual_3_sensors_high.mat'); % SSI modes
K_SSI_mode_sens = dataSSImode_sens.knew';

dataSSIFreqmode_sens = load('.\Sensitivity_method\data_updated_par_sens_3_sensors\Eigenvalue_Mode_shape_residual_3_sensors_high.mat'); % SSI modes and freq
K_SSI_freq_mode_sens = dataSSIFreqmode_sens.knew';


% PLotting the stiffness
y = [K_SSI_freq;K_SSI_freq_sens;K_SSI_mode;K_SSI_mode_mac;K_SSI_mode_sens;K_SSI_freq_mode;K_SSI_freq_mode_sens;...
    K_SSI_freqmodeEIL;Globalstiff];
figure
hold on
b = bar(y,'stacked');
for i = 1:size(y,2)
    xtips = b(i).XEndPoints;
    ytips = b(i).YEndPoints;
    labels = string(b(i).YData);
    text(xtips,ytips./1.1,labels,'HorizontalAlignment','center','VerticalAlignment','middle')
end
hold off
legend('k1','k2','k3','k4','k5')
title('Values of stiffness for different OMA methods for cost functions and sensitivity method','FontSize',20)
xlabel('Method')
xticks(xtips)
xticklabels({'(freq)','freq sens',' (mode)',' (mode mac)','mode sens',' (freq+mode)','freq+mode sens',...
    ' EIL (freq+mode)','Global stiffness','FontSize',14})
ylabel('Stiffness [N/m]','FontSize',14)

% 3D plot of the stiffness
figure
bar3(y)
title('Values of stiffness for different OMA methods for cost functions and sensitivity method','FontSize',20)
%xlabel('ks')
%ylabel('Method')
xticks([1,2,3,4,5])
xticklabels({'k1','k2','k3','k4','k5','FontSize',14})
yticks(xtips)
yticklabels({'(freq)','freq sens',' (mode)',' (mode mac)','mode sens',' (freq+mode)','freq+mode sens',...
    ' EIL (freq+mode)','Global stiffness','FontSize',14})
zlabel('Stiffness [N/m]','FontSize',14)


