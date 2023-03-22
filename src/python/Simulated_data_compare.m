
% parameters
clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'),genpath('python'),genpath('npy-matlab-master'))
% Loading modal parameters from OMA 
SSIFreq = readNPY('omega.npy');
SSIomega = SSIFreq * 2 * pi;
SSIphi = readNPY('phi.npy');
SSIdamp = readNPY('damp.npy');

SSIdatFreq = readNPY('SSIdatomega.npy');
SSIdatomega = SSIFreq * 2 * pi;
SSIdatphi = readNPY('SSIdatphi.npy');
SSIdatdamp = readNPY('SSIdatdamp.npy');

FDDFreq = readNPY('FDDomega.npy');
FDDomega = FDDFreq * 2 * pi;
FDDphi = readNPY('FDDphi.npy');
FDDdamp = readNPY('FDDdamp.npy');

% numerical
filename = load('modelprop_jan.mat'); % omegas from numericla model
fn = filename.fn;
U = filename.U;

% simulated damping
sim = load('data_sim_newmark_jan_damp.mat');
simdamp = sim.zetas;

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

% damping
modes = {'Mode 1';'Mode 2';'Mode 3';'Mode 4';'Mode 5'};
damp = table(simdamp,SSIdamp,SSIdatdamp,FDDdamp,...
    'RowNames',modes);
disp('Damping from each method compared to simulated :')
disp(damp)

prompt = "Which OMA is used (1=SSIcov, 2 = SSIdat, 3=FDD)? ";
MODE = input(prompt);
if MODE == 1
    OMAfreq = SSIFreq;
    OMAphi = SSIphi;
elseif MODE == 2
    OMAfreq = SSIdatFreq;
    OMAphi = SSIdatphi;
elseif MODE == 3
    OMAfreq = FDDFreq;
    OMAphi = FDDphi;
else
    disp('ADAM! WHAT DOES HØVL MEAN???')
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
for i=1:length(fn)
    subplot(1,length(fn),i)
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
    title(['f = ' num2str(OMAfreq(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
    xline(0.0,'--')
    xlim([-1.1,1.1])
    ylim([0,x(end)])
    if i==1 && MODE==1
        legend('Numerical','SSIcov','Location','northwest')
    elseif i==1 && MODE==2
        legend('Numerical','SSIDat','Location','northwest')
    elseif i==1
        legend('Numerical','FDD','Location','northwest')
    end
end
if MODE == 1
    sgtitle('Numerical mode shapes, calibrated by SSIcov','FontSize',20)
elseif MODE == 2
    sgtitle('Numerical mode shapes, calibrated by SSIdat','FontSize',20)
elseif MODE == 3
    sgtitle('Numerical mode shapes, calibrated by FDD','FontSize',20)
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);

% Frequency accaruancy
disp(strcat('Frequency accuracy,1 : ',num2str(min(OMAfreq(1),fn(1))/max(OMAfreq(1),fn(1))*100),'%'));
disp(strcat('Frequency accuracy,2 : ',num2str(min(OMAfreq(2),fn(2))/max(OMAfreq(2),fn(2))*100),'%'));
disp(strcat('Frequency accuracy,3 : ',num2str(min(OMAfreq(3),fn(3))/max(OMAfreq(3),fn(3))*100),'%'));
disp(strcat('Frequency accuracy,4 : ',num2str(min(OMAfreq(4),fn(4))/max(OMAfreq(4),fn(4))*100),'%'));
disp(strcat('Frequency accuracy,5 : ',num2str(min(OMAfreq(5),fn(5))/max(OMAfreq(5),fn(5))*100),'%'));
disp(strcat('Mean frequency accuracy : ',num2str(mean(min(OMAfreq,fn)./max(OMAfreq,fn)*100)),'%'));

% CrossMAC plot of mode shapes
mac=crossMAC(U,OMAphi,MODE,[OMAfreq,fn]);
dmac = diag(mac);
if MODE==1
    disp('Modal Assurance Criterion between Numerical modeshapes and SSIcov  : ')
    disp(strcat(num2str(mac)));
elseif MODE==2
    disp('Modal Assurance Criterion between Numerical modeshapes and SSIdat  : ')
    disp(strcat(num2str(mac)));
elseif MODE == 3
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

