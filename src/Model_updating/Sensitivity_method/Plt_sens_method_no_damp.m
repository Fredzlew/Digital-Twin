%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting the model update no damp %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath(genpath('..\..\data'),genpath('.\data_updated_par_sens'),genpath('..\Costfunctions\functions'))


% Loading stiffness for all OMA sens and numerical model
filename = load('..\..\data\modelprop.mat'); % omegas from numericla model
omegas = filename.fn * 2 * pi;
Globalstiff = filename.k;
H = filename.H;
Km = filename.K;

dataSSIFreq = load('.\data_updated_par_sens\Eigenvalue_residual_no_damp.mat'); % SSI FREQ
K_SSI_freq = dataSSIFreq.knew';

dataSSImode = load('.\data_updated_par_sens\Mode_shape_residual_no_damp.mat'); % SSI modes
K_SSI_mode = dataSSImode.knew';

dataSSIFreqmode = load('.\data_updated_par_sens\Eigenvalue_Mode_shape_residual_no_damp.mat'); % SSI modes and freq
K_SSI_freq_mode = dataSSIFreqmode.knew';




% PLotting the stiffness
y = [K_SSI_freq;K_SSI_mode;K_SSI_freq_mode;...
    Globalstiff];
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
title('Values of stiffness for different OMA methods and sensitivity method','FontSize',20)
xlabel('Method')
xticks(xtips)
xticklabels({'SSI (freq)','SSI (mode)','SSI (freq+mode)',...
    'Global stiffness','FontSize',14})
ylabel('Stiffness [N/m]','FontSize',14)

% 3D plot of the stiffness
figure
bar3(y)
title('Values of stiffness for different OMA methods and sensitivity method','FontSize',20)
%xlabel('ks')
%ylabel('Method')
xticks([1,2,3,4,5])
xticklabels({'k1','k2','k3','k4','k5','FontSize',14})
yticks(xtips)
yticklabels({'freq','mode','freq+mode',...
    'Global stiffness','FontSize',14})
zlabel('Stiffness [N/m]','FontSize',14)

% What we want to plot
% measuered natural frequencies from OMA SSI-cov [Hz]
OMAfreq = readNPY('..\..\data\experimental_data\Modal_par\SSIfreq_no_damp.npy');

% Load mode shapes
OMAphi = readNPY('..\..\data\experimental_data\Modal_par\SSImodes_no_damp.npy');


promptt = "Which do you want to plot? (1=SSI (freq), 2=SSI (mode) " + ...
    " 3=SSI (freq+mode)? ";
xx = input(promptt);
if xx == 1
    data = load('.\data_updated_par_sens\Eigenvalue_residual_no_damp.mat'); % SSI FREQ
    fn = data.fn;
    K = data.Knew;
    U = data.U;
elseif xx == 2
    data = load('.\data_updated_par_sens\Mode_shape_residual_no_damp.mat'); % SSI modes
    fn = data.fn;
    K = data.Knew;
    U = data.U;
elseif xx == 3
    data = load('.\data_updated_par_sens\Eigenvalue_Mode_shape_residual_no_damp.mat'); % SSI EIL  FREQ and modes
    fn = data.fn;
    K = data.Knew;
    U = data.U;
end

disp('----------------------------------------------------------------------')

% plotting the mode shapes
x = [0, H];
phi = [zeros(1,length(U)); U];
fig = figure;
fig.Position=[100 100 1600 700];

for i = 1:5
    if phi(2,i)*OMAphi(1,i) < 0 
        phi(:,i) = phi(:,i)*-1;
    end
end

for i=1:length(omegas)
    subplot(1,length(omegas),i)
    hold on
    plot(phi(:,i),x,'-m')
    if phi(2,i)*OMAphi(1,i) < 0 % Swap sign on mode shape
        plot([0  ;-OMAphi(:,i)],x,'go-.');
        plot(-OMAphi(1:end,i),x(2:end),'g.','markersize',30)
    else
        plot([0  ;OMAphi(:,i)],x,'go-.');
        plot(OMAphi(1:end,i),x(2:end),'g.','markersize',30)
    end
    plot(phi(2:end,i),x(2:end),'b.','markersize',30)
    title(['f = ' num2str(fn(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
    xline(0.0,'--')
    xlim([-1.1,1.1])
    ylim([0,x(end)])
    if i==1 
        legend('Numerical','SSI','Location','northwest')
    end
end

sgtitle('Numerical mode shapes, calibrated by SSI','FontSize',20)

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);

disp(strcat('Frequency accuracy,1 : ',num2str(min(OMAfreq(1),fn(1))/max(OMAfreq(1),fn(1))*100),'%'));
disp(strcat('Frequency accuracy,2 : ',num2str(min(OMAfreq(2),fn(2))/max(OMAfreq(2),fn(2))*100),'%'));
disp(strcat('Frequency accuracy,3 : ',num2str(min(OMAfreq(3),fn(3))/max(OMAfreq(3),fn(3))*100),'%'));
disp(strcat('Frequency accuracy,4 : ',num2str(min(OMAfreq(4),fn(4))/max(OMAfreq(4),fn(4))*100),'%'));
disp(strcat('Frequency accuracy,5 : ',num2str(min(OMAfreq(5),fn(5))/max(OMAfreq(5),fn(5))*100),'%'));
disp(strcat('Mean frequency accuracy : ',num2str(mean(min(OMAfreq,fn)./max(OMAfreq,fn)*100)),'%'));


% Display stiffness matrices and changes in these
disp('----------------------------------------------------------------------')
disp('Original stiffness matrix [N/m]  : ')
disp(strcat(num2str(Km)));
disp('Calibrated stiffness matrix [N/m]  : ')
disp(strcat(num2str(K)));
disp('Changes in stiffness matrix [N/m] : ')
disp(strcat(num2str(abs(abs(K)-abs(Km)))));
disp('Total change in stiffness matrix [N/m] : ')
disp(strcat(num2str(sum(sum(abs(abs(K)-abs(Km)))))));
disp('----------------------------------------------------------------------')

% CrossMAC plot of mode shapes
mac=crossMAC(U,OMAphi);
dmac = diag(mac);

disp('Modal Assurance Criterion between Numerical modeshapes and SSI  : ')
disp(strcat(num2str(mac)));
disp('----------------------------------------------------------------------')
disp(strcat('Mode shape accuracy (MAC),1 : ',num2str(dmac(1)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),2 : ',num2str(dmac(2)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),3 : ',num2str(dmac(3)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),4 : ',num2str(dmac(4)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),5 : ',num2str(dmac(5)*100),'%'));
disp(strcat('Mean mode shape accuracy (MAC): ',num2str(mean(dmac)*100),'%'));
disp('----------------------------------------------------------------------')

% % Plot MAC
figure
barMAC = bar3(mac);
for k = 1:length(barMAC)   
    zdata = barMAC(k).ZData;   
    barMAC(k).CData = zdata;
    barMAC(k).FaceColor = 'interp';
end
colormap(jet);
colorbar
title('MAC - Numerical compared to SSI no_damp damp')
xlabel('FE-model frequencies [Hz]')
ylabel('OMA frequencies [Hz]')
xticks([1,2,3,4,5])
xticklabels(string(fn'))
yticks([1,2,3,4,5])
yticklabels(string(OMAfreq'))
box on

%% Download plot
% PLotting the stiffness
if xx == 1
k1 = y(:,1);
k2 = y(:,2);
k3 = y(:,3);
k4 = y(:,4);
k5 = y(:,5);
%T_stiff = array2table([num2cell(y)]);
T_modeshapes = array2table([num2cell(x'),num2cell(phi),num2cell([0;fn])]);
T_mac = array2table([num2cell(mac)]);

stiff= {'Freq';'Mode';'Freq mode';'Global stiffness'};
T_stiff = table(k1,k2,k3,k4,k5, 'RowNames',stiff);

%T_stiff.Properties.VariableNames(1:6) = {'','k1','k2','k3','k4','k5'};
T_modeshapes.Properties.VariableNames(1:7) = {'height','numphi1','numphi2','numphi3','numphi4','numphi5','numfreq'};
T_mac.Properties.VariableNames(1:5) = {'mac1','mac2','mac3','mac4','mac5'};


writetable(T_stiff,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap9_sens_stiff_nodamp.csv','Delimiter',';')
writetable(T_modeshapes,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap9_sens_modeshapes1_nodamp.csv','Delimiter',';')
writetable(T_mac,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap9_sens_mac1_nodamp.csv','Delimiter',';')

elseif xx == 2
T_modeshapes = array2table([num2cell(x'),num2cell(phi),num2cell([0;fn])]);
T_mac = array2table([num2cell(mac)]);

T_modeshapes.Properties.VariableNames(1:7) = {'height','numphi1','numphi2','numphi3','numphi4','numphi5','numfreq'};
T_mac.Properties.VariableNames(1:5) = {'mac1','mac2','mac3','mac4','mac5'};

writetable(T_modeshapes,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap9_sens_modeshapes2_nodamp.csv','Delimiter',';')
writetable(T_mac,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap9_sens_mac2_nodamp.csv','Delimiter',';')
elseif xx == 3
T_modeshapes = array2table([num2cell(x'),num2cell(phi),num2cell([0;fn])]);
T_mac = array2table([num2cell(mac)]);

T_modeshapes.Properties.VariableNames(1:7) = {'height','numphi1','numphi2','numphi3','numphi4','numphi5','numfreq'};
T_mac.Properties.VariableNames(1:5) = {'mac1','mac2','mac3','mac4','mac5'};

writetable(T_modeshapes,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap9_sens_modeshapes3_nodamp.csv','Delimiter',';')
writetable(T_mac,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap9_sens_mac3_nodamp.csv','Delimiter',';')
end