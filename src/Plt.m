%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting the modal update %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))

% Loading stiffness for all OMA costfunction and numerical model
filename = load('modelprop.mat'); % omegas from numericla model
omegas = filename.fn * 2 * pi;

dataSSIFreq = load('costfunupdateSSIfreq.mat'); % SSI FREQ
K_SSI_freq = diag(dataSSIFreq.K)';
H = dataSSIFreq.H;
k2 = dataSSIFreq.k2;
U = dataSSIFreq.U;

dataERAFreq = load('costfunupdateERAfreq.mat'); % ERA FREQ
K_ERA_freq = diag(dataERAFreq.K)';

dataFDDFreq = load('costfunupdateFDDfreq.mat'); % FDD FREQ
K_FDD_freq = diag(dataERAFreq.K)';

dataSSImode = load('costfunupdateSSImode.mat'); % SSI modes
K_SSI_mode = diag(dataSSImode.K)';

dataERAmode = load('costfunupdateERAmode.mat'); % ERA modes
K_ERA_mode = diag(dataERAmode.K)';

dataFDDmode = load('costfunupdateFDDmode.mat'); % FDD modes
K_FDD_mode = diag(dataERAmode.K)';

dataSSIFreqmode = load('costfunupdateSSIfreqmode.mat'); % SSI FREQ and modes
K_SSI_freq_mode = diag(dataSSIFreqmode.K)';

dataERAFreqmode = load('costfunupdateERAfreqmode.mat'); % ERA FREQ and modes
K_ERA_freq_mode = diag(dataERAFreqmode.K)';

dataFDDFreqmode = load('costfunupdateFDDfreqmode.mat'); % FDD FREQ and modes
K_FDD_freq_mode = diag(dataERAFreqmode.K)';

dataSSIfreqmodeEIL = load('costfunupdateSSIfreqmodeEIL.mat'); % SSI modes
K_SSI_freqmodeEIL = diag(dataSSIfreqmodeEIL.K)';

dataERAfreqmodeEIL = load('costfunupdateERAfreqmodeEIL.mat'); % ERA modes
K_ERA_freqmodeEIL = diag(dataERAfreqmodeEIL.K)';

dataFDDfreqmodeEIL = load('costfunupdateFDDfreqmodeEIL.mat'); % FDD modes
K_FDD_freqmodeEIL = diag(dataFDDfreqmodeEIL.K)';

dataSSIfreqmodeEILJAN = load('costfunupdateSSIfreqmodeEILJAN.mat'); % SSI modes
K_SSI_freqmodeEILJAN = diag(dataSSIfreqmodeEIL.K)';

dataERAfreqmodeEILJAN = load('costfunupdateERAfreqmodeEILJAN.mat'); % ERA modes
K_ERA_freqmodeEILJAN = diag(dataERAfreqmodeEILJAN.K)';

dataFDDfreqmodeEILJAN = load('costfunupdateFDDfreqmodeEILJAN.mat'); % FDD modes
K_FDD_freqmodeEILJAN = diag(dataFDDfreqmodeEILJAN.K)';



% PLotting the stiffness
y = [K_SSI_freq;K_ERA_freq;K_FDD_freq;K_SSI_mode;K_ERA_mode;K_FDD_mode;K_SSI_freq_mode;K_ERA_freq_mode;K_FDD_freq_mode;...
    K_SSI_freqmodeEIL;K_ERA_freqmodeEIL;K_FDD_freqmodeEIL;K_SSI_freqmodeEILJAN;K_ERA_freqmodeEILJAN;K_FDD_freqmodeEILJAN;k2];
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
title('Values of stiffness for different OMA methods and cost functions')
xlabel('Method')
xticks(xtips)
xticklabels({'SSI (freq)','ERA (freq)','FDD (freq)','SSI (mode)','ERA (mode)','FDD (mode)','SSI (freq+mode)','ERA (freq+mode)','FDD (freq+mode)',...
    'SSI EIL (freq+mode)','ERA EIL (freq+mode)','FDD EIL (freq+mode)','SSI EIL JAN(freq+mode)','ERA EIL JAN(freq+mode)','FDD EIL JAN(freq+mode)','Geometric Stiffness'})
ylabel('Stiffness [N/m]')

% 3D plot of the stiffness
figure
bar3(y)
title('Values of stiffness for different OMA methods and cost functions')
%xlabel('ks')
%ylabel('Method')
xticks([1,2,3,4,5])
xticklabels({'k1','k2','k3','k4','k5'})
yticks(xtips)
yticklabels({'SSI (freq)','ERA (freq)','FDD (freq)','SSI (mode)','ERA (mode)','FDD (mode)','SSI (freq+mode)','ERA (freq+mode)','FDD (freq+mode)',...
    'SSI EIL (freq+mode)','ERA EIL (freq+mode)','FDD EIL (freq+mode)','SSI EIL JAN(freq+mode)','ERA EIL JAN(freq+mode)','FDD EIL JAN(freq+mode)','Geometric Stiffness'})
zlabel('Stiffness [N/m]')

% What we want to plot
% Define what OMA method is used 
% MODE = 2; % 1=SSI, 2=ERA, 3=FDD
prompt = "Which OMA is used (1=SSI, 2=ERA, 3=FDD)? ";
MODE = input(prompt);
promptt = "Which do you want to plot? (1=SSI (freq), 2=ERA (freq), 3=FDD (freq), 4=SSI (mode), 5=ERA (mode), 6=FDD (mode)," + ...
    " 7=SSI (freq+mode), 8=ERA (freq+mode), 9=FDD (freq+mode), 10=SSI EIL (freq+mode), 11=ERA EIL (freq+mode), 12=FDD EIL (freq+mode)," + ...
    " 13=SSI EIL JAN(freq+mode), 14=ERA EIL JAN(freq+mode), 15=FDD EIL JAN(freq+mode))? ";
x = input(promptt);
if x == 1
    data = load('costfunupdateSSIfreq.mat'); % SSI FREQ
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 2
    data = load('costfunupdateERAfreq.mat'); % ERA FREQ
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 3
    data = load('costfunupdateFDDfreq.mat'); % FDD FREQ
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 4
    data = load('costfunupdateSSImode.mat'); % SSI modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 5
    data = load('costfunupdateERAmode.mat'); % ERA modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 6
    data = load('costfunupdateFDDmode.mat'); % FDD modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 7
    data = load('costfunupdateSSIfreqmode.mat'); % SSI FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 8
    data = load('costfunupdateERAfreqmode.mat'); % ERA FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 9
    data = load('costfunupdateFDDfreqmode.mat'); % FDD FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 10 
    data = load('costfunupdateSSIfreqmodeEIL.mat'); % SSI EIL FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 11
    data = load('costfunupdateERAfreqmodeEIL.mat'); % ERA EIL FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 12
    data = load('costfunupdateFDDfreqmodeEIL.mat'); % FDD EIL FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 13
    data = load('costfunupdateSSIfreqmodeEILJAN.mat'); % SSI EIL JAN FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 14
    data = load('costfunupdateERAfreqmodeEILJAN.mat'); % ERA EIL JAN FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
elseif x == 15
    data = load('costfunupdateFDDfreqmodeEILJAN.mat'); % FDD EIL JAN FREQ and modes
    OMAphi = data.OMAphi;
    OMAfreq = data.OMAfreq;
    fn = data.fn;
    K = data.K;
    Km = data.Km;
end

% plotting the mode shapes
x = [0, H];
phi = [zeros(1,length(U)); U];
fig = figure;
fig.Position=[100 100 1600 700];
% if MODE==1
%     OMAphi = SSIphi;
% elseif MODE==2
%     OMAphi = ERAphi;
% else
%     OMAphi = FDDphi;
% end
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
    if i==1 && MODE==1
        legend('Numerical','SSI','Location','northwest')
    elseif i==1 && MODE==2
        legend('Numerical','ERA','Location','northwest')
    elseif i==1
        legend('Numerical','FDD','Location','northwest')
    end
end
if MODE==1
    sgtitle('Numerical mode shapes, calibrated by SSI','FontSize',20)
elseif MODE==2
    sgtitle('Numerical mode shapes, calibrated by ERA','FontSize',20)
else
    sgtitle('Numerical mode shapes, calibrated by FDD','FontSize',20)
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);

% Display accuracy of natural frequencies
% if MODE==1
%     OMAfreq=SSIFreq;
%     OMAphi=SSIphi;
% elseif MODE==2
%     OMAfreq=ERAFreq;
%     OMAphi=ERAphi;
% else
%     OMAfreq=FDDFreq;
%     OMAphi=FDDphi;
% end


disp(strcat('Frequency accuracy,1 : ',num2str(min(OMAfreq(1),fn(1))/max(OMAfreq(1),fn(1))*100),'%'));
disp(strcat('Frequency accuracy,2 : ',num2str(min(OMAfreq(2),fn(2))/max(OMAfreq(2),fn(2))*100),'%'));
disp(strcat('Frequency accuracy,3 : ',num2str(min(OMAfreq(3),fn(3))/max(OMAfreq(3),fn(3))*100),'%'));
disp(strcat('Frequency accuracy,4 : ',num2str(min(OMAfreq(4),fn(4))/max(OMAfreq(4),fn(4))*100),'%'));
disp(strcat('Frequency accuracy,5 : ',num2str(min(OMAfreq(5),fn(5))/max(OMAfreq(5),fn(5))*100),'%'));
disp(strcat('Mean frequency accuracy : ',num2str(mean(min(OMAfreq,fn)./max(OMAfreq,fn)*100)),'%'));

% Display accuracy of mode shapes
[MSacc,TOTacc]=modeshapeacc(OMAphi,U);
disp('----------------------------------------------------------------------')
disp(strcat('Mode shape accuracy,1 : ',num2str(MSacc(1)*100),'%'));
disp(strcat('Mode shape accuracy,2 : ',num2str(MSacc(2)*100),'%'));
disp(strcat('Mode shape accuracy,3 : ',num2str(MSacc(3)*100),'%'));
disp(strcat('Mode shape accuracy,4 : ',num2str(MSacc(4)*100),'%'));
disp(strcat('Mode shape accuracy,5 : ',num2str(MSacc(5)*100),'%'));
disp(strcat('Mean mode shape accuracy : ',num2str(TOTacc*100),'%'));


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
disp(strcat('Total mean accuracy (mode+freq) : ',num2str(mean([TOTacc*100,mean(min(OMAfreq,fn)./max(OMAfreq,fn)*100)])),'%'));
disp('----------------------------------------------------------------------')

% CrossMAC plot of mode shapes
mac=crossMAC(U,OMAphi,MODE,[OMAfreq,fn]);
dmac = diag(mac);
if MODE==1
    disp('Modal Assurance Criterion between Numerical modeshapes and SSI  : ')
    disp(strcat(num2str(mac)));
elseif MODE==2
    disp('Modal Assurance Criterion between Numerical modeshapes and ERA  : ')
    disp(strcat(num2str(mac)));
else
    disp('Modal Assurance Criterion between Numerical modeshapes and FDD  :' )
    disp(strcat(num2str(mac)));
end
disp('----------------------------------------------------------------------')
disp(strcat('Mode shape accuracy (MAC),1 : ',num2str(dmac(1)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),2 : ',num2str(dmac(2)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),3 : ',num2str(dmac(3)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),4 : ',num2str(dmac(4)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),5 : ',num2str(dmac(5)*100),'%'));
disp(strcat('Mean mode shape accuracy (MAC): ',num2str(mean(dmac)*100),'%'));
disp('----------------------------------------------------------------------')