%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting the model update high %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath(genpath('..\..\data'),genpath('.\data_updated_par'),genpath('.\functions'))

% Loading stiffness for all OMA costfunction and numerical model
filename = load('..\..\data\modelprop.mat'); % omegas from numericla model
omegas = filename.fn * 2 * pi;
L_fe = filename.L;
EI_fe = filename.EI;
H = filename.H;
Globalstiff = filename.k;

dataSSIFreq = load('.\data_updated_par\SSIfreq_high.mat'); % SSI FREQ
K_SSI_freq = dataSSIFreq.k;

dataSSImode = load('.\data_updated_par\SSImode_high.mat'); % SSI modes
K_SSI_mode = dataSSImode.k;

dataSSImodemac = load('.\data_updated_par\SSImode_mac_high.mat'); % SSI modes
K_SSI_mode_mac = dataSSImodemac.k;


dataSSIFreqmode = load('.\data_updated_par\SSIfreqmode_high.mat'); % SSI FREQ and modes
K_SSI_freq_mode = dataSSIFreqmode.k;

dataSSIfreqmodeEIL = load('.\data_updated_par\SSIfreqmodeEIL_high.mat'); % SSI modes
K_SSI_freqmodeEIL = dataSSIfreqmodeEIL.k;
L = dataSSIfreqmodeEIL.L;
EI = dataSSIfreqmodeEIL.EI;

% PLotting the stiffness
y = [K_SSI_freq;K_SSI_mode;K_SSI_mode_mac;K_SSI_freq_mode;...
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
title('Values of stiffness for different OMA methods and cost functions','FontSize',20)
xlabel('Method')
xticks(xtips)
xticklabels({'SSI (freq)','SSI (mode)','SSI (mode mac)','SSI (freq+mode)',...
    'SSI EIL (freq+mode)','Global stiffness','FontSize',14})
ylabel('Stiffness [N/m]','FontSize',14)

% 3D plot of the stiffness
figure
bar3(y)
title('Values of stiffness for different OMA methods and cost functions','FontSize',20)
%xlabel('ks')
%ylabel('Method')
xticks([1,2,3,4,5])
xticklabels({'k1','k2','k3','k4','k5','FontSize',14})
yticks(xtips)
yticklabels({'SSI (freq)','SSI (mode)','SSI (mode mac)','SSI (freq+mode)',...
    'SSI EIL (freq+mode)','Global stiffness','FontSize',14})
zlabel('Stiffness [N/m]','FontSize',14)

% What we want to plot
% Define what OMA method is used 
% MODE = 2; % 1=SSI, 2=ERA, 3=FDD

promptt = "Which do you want to plot? (1=SSI (freq), 2=SSI (mode), 3=SSI (mode mac) " + ...
    " 4=SSI (freq+mode), " + ...
    " 5=SSI EIL (freq+mode)? ";
xx = input(promptt);
if xx == 1
    data = load('.\data_updated_par\SSIfreq_high.mat'); % SSI FREQ
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
    U = data.U;
elseif xx == 2
    data = load('.\data_updated_par\SSImode_high.mat'); % SSI modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
    U = data.U;
elseif xx == 3
    data = load('.\data_updated_par\SSImode_mac_high.mat'); % SSI EIL  FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
    U = data.U;
elseif xx == 4
    data = load('.\data_updated_par\SSIfreqmode_high.mat'); % SSI FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
    U = data.U;
elseif xx == 5
    data = load('.\data_updated_par\SSIfreqmodeEIL_high.mat'); % SSI EIL  FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
    U = data.U;
    L = data.L;
    EI = data.EI;
end
if xx == 5
    disp(strcat('Height : ',num2str(L_fe),'m'));
    disp(strcat('EI  : ',num2str(EI_fe),'Nm^2'));
    disp(strcat('Height after costfunction : ',num2str(L),'m'));
    disp(strcat('EI after costfunction : ',num2str(EI),'Nm^2'));
    disp(strcat('Height difference : ',num2str(L-L_fe),'m'));
    disp(strcat('EI difference : ',num2str(EI-EI_fe),'Nm^2'));
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
title('MAC - Numerical compared to SSI high damp')
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
interval = [1;2;3;4;5;6];
k1 = y(:,1);
k2 = y(:,2);
k3 = y(:,3);
k4 = y(:,4);
k5 = y(:,5);
%T_stiff = array2table([num2cell(y)]);
T_modeshapes = array2table([num2cell(x'),num2cell(phi),num2cell([0;fn])]);
T_mac = array2table([num2cell(mac)]);

stiff= {'Freq';'Mode';'Mode mac';'Freq + mode';'Fred + mode EI L';'Global stiffness'};
T_stiff = table(interval,k1,k2,k3,k4,k5, 'RowNames',stiff);

%T_stiff.Properties.VariableNames(1:6) = {'','k1','k2','k3','k4','k5'};
T_modeshapes.Properties.VariableNames(1:7) = {'height','numphi1','numphi2','numphi3','numphi4','numphi5','numfreq'};
T_mac.Properties.VariableNames(1:5) = {'mac1','mac2','mac3','mac4','mac5'};
writetable(T_stiff,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_stiff_highdamp.csv','Delimiter',';')
writetable(T_modeshapes,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_modeshapes1_highdamp.csv','Delimiter',';')
writetable(T_mac,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_mac1_highdamp.csv','Delimiter',';')
elseif xx == 2
T_modeshapes = array2table([num2cell(x'),num2cell(phi),num2cell([0;fn])]);
T_mac = array2table([num2cell(mac)]);

T_modeshapes.Properties.VariableNames(1:7) = {'height','numphi1','numphi2','numphi3','numphi4','numphi5','numfreq'};
T_mac.Properties.VariableNames(1:5) = {'mac1','mac2','mac3','mac4','mac5'};

writetable(T_modeshapes,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_modeshapes2_highdamp.csv','Delimiter',';')
writetable(T_mac,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_mac2_highdamp.csv','Delimiter',';')
elseif xx == 3
T_modeshapes = array2table([num2cell(x'),num2cell(phi),num2cell([0;fn])]);
T_mac = array2table([num2cell(mac)]);

T_modeshapes.Properties.VariableNames(1:7) = {'height','numphi1','numphi2','numphi3','numphi4','numphi5','numfreq'};
T_mac.Properties.VariableNames(1:5) = {'mac1','mac2','mac3','mac4','mac5'};

writetable(T_modeshapes,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_modeshapes3_highdamp.csv','Delimiter',';')
writetable(T_mac,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_mac3_highdamp.csv','Delimiter',';')
elseif xx == 4
T_modeshapes = array2table([num2cell(x'),num2cell(phi),num2cell([0;fn])]);
T_mac = array2table([num2cell(mac)]);

T_modeshapes.Properties.VariableNames(1:7) = {'height','numphi1','numphi2','numphi3','numphi4','numphi5','numfreq'};
T_mac.Properties.VariableNames(1:5) = {'mac1','mac2','mac3','mac4','mac5'};

writetable(T_modeshapes,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_modeshapes4_highdamp.csv','Delimiter',';')
writetable(T_mac,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_mac4_highdamp.csv','Delimiter',';')
elseif xx == 5
T_modeshapes = array2table([num2cell(x'),num2cell(phi),num2cell([0;fn])]);
T_mac = array2table([num2cell(mac)]);

T_modeshapes.Properties.VariableNames(1:7) = {'height','numphi1','numphi2','numphi3','numphi4','numphi5','numfreq'};
T_mac.Properties.VariableNames(1:5) = {'mac1','mac2','mac3','mac4','mac5'};

writetable(T_modeshapes,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_modeshapes5_highdamp.csv','Delimiter',';')
writetable(T_mac,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap8_costfunction_mac5_highdamp.csv','Delimiter',';')
end
