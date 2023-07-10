
% parameters
clc; clear; close all;
addpath(genpath('..\..\data'),genpath('..\..\npy-matlab-master'))
% Loading modal parameters from OMA simulation
SSIphi = readNPY('..\..\data\simulated_data\Modal_par\SSIcovmodes.npy');
SSIFreq = readNPY('..\..\data\simulated_data\Modal_par\SSIcovfreq.npy');
SSIomega = SSIFreq * 2 * pi;
SSIdamp = readNPY('..\..\data\simulated_data\Modal_par\SSIcovdamp.npy');

SSIdatFreq = readNPY('..\..\data\simulated_data\Modal_par\SSIdatfreq.npy');
SSIdatomega = SSIFreq * 2 * pi;
SSIdatphi = readNPY('..\..\data\simulated_data\Modal_par\SSIdatmodes.npy');
SSIdatdamp = readNPY('..\..\data\simulated_data\Modal_par\SSIdatdamp.npy');

FDDFreq = readNPY('..\..\data\simulated_data\Modal_par\FDDfreq.npy');
FDDomega = FDDFreq * 2 * pi;
FDDphi = readNPY('..\..\data\simulated_data\Modal_par\FDDmodes.npy');
FDDdamp = readNPY('..\..\data\simulated_data\Modal_par\FDDdamp.npy');

% number of modes
nm = 5;

% numerical
filename = load('..\..\data\modelprop.mat'); % omegas from numericla model
fn = filename.fn;
U = filename.U;

% simulated damping
for i = 1:1
sim = load([num2str(i) '_sim_data.mat']);
simdamp1(:,i) = sim.zetas;
end
simdamp = (simdamp1);


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

SSIcovfreq = zeros(nm,length(SSIFreq)/nm);
FDDfreq = zeros(nm,length(FDDFreq)/nm);
SSIcovdamp = zeros(nm,length(SSIdamp)/nm);
FDDdamping = zeros(nm,length(FDDdamp)/nm);
% Setting it up in the right order 
for i = 1:length(SSIFreq)/nm
    for j = 1:5
    SSIcovfreq(j,i) = SSIFreq(i+(j-1)+(i-1)*4);
    FDDfreq(j,i) = FDDFreq(i+(j-1)+(i-1)*4);
    SSIcovdamp(j,i) = SSIdamp(i+(j-1)+(i-1)*4);
    FDDdamping(j,i) = FDDdamp(i+(j-1)+(i-1)*4);
    end
end

% Mean value and standard deviation
SSIcovmu_freq = zeros(nm,1);
FDDmu_freq = zeros(nm,1);
SSIcovmu_damp = zeros(nm,1);
FDDmu_damp = zeros(nm,1);

SSIcovsd_freq = zeros(nm,1);
FDDsd_freq = zeros(nm,1);
SSIcovsd_Damp = zeros(nm,1);
FDDsd_Damp = zeros(nm,1);
for i = 1:nm
    % mean
    SSIcovmu_freq(i) = mean(SSIcovfreq(i,:));
    FDDmu_freq(i) = mean(FDDfreq(i,:));
    SSIcovmu_damp(i) = mean(SSIcovdamp(i,:));
    FDDmu_damp(i) = mean(FDDdamping(i,:));
    % standard deviation
    SSIcovsd_freq(i) = std(SSIcovfreq(i,:));
    FDDsd_freq(i) = std(FDDfreq(i,:));
    SSIcovsd_Damp(i) = std(SSIcovdamp(i,:));
    FDDsd_Damp(i) = std(FDDdamping(i,:));
end

% mean and standard deviation for modes
SSIcovmu_modes = zeros(nm,nm);
FDDmu_modes = zeros(nm,nm);
SSIcovsd_modes = zeros(nm,nm);
FDDsd_modes = zeros(nm,nm);
for i = 1:nm
    for j = 1:nm
        % mean
        SSIcovmu_modes(i,j) = mean(SSIphi(:,j,i));
        FDDmu_modes(i,j) = mean(FDDphi(:,j,i));
        % standard deviation
        SSIcovsd_modes(i,j) = std(SSIphi(:,j,i));
        FDDsd_modes(i,j) = std(FDDphi(:,j,i));
    end
end






% damping
modes = {'Mode 1';'Mode 2';'Mode 3';'Mode 4';'Mode 5'};
damp = table(simdamp,SSIcovmu_damp,SSIdatdamp,FDDmu_damp,...
    'RowNames',modes);
disp('Mean damping from each method over 1000 simulations :')
disp(damp)

% damping
modes = {'Mode 1';'Mode 2';'Mode 3';'Mode 4';'Mode 5'};
damp = table(SSIcovsd_Damp,FDDsd_Damp,...
    'RowNames',modes);
disp('standard deviation damping from each method over 1000 simulations :')
disp(damp)

% frekvens
modes = {'Mode 1';'Mode 2';'Mode 3';'Mode 4';'Mode 5'};
damp = table(fn,SSIcovmu_freq,SSIdatFreq,FDDmu_freq,...
    'RowNames',modes);
disp('Mean frequenzy from each method over 1000 simulations :')
disp(damp)

% frekvens
modes = {'Mode 1';'Mode 2';'Mode 3';'Mode 4';'Mode 5'};
damp = table(SSIcovsd_freq,FDDsd_freq,...
    'RowNames',modes);
disp('standard deviation frequenzy from each method over 1000 simulations :')
disp(damp)


prompt = "Which OMA is used (1=SSIcov, 2 = SSIdat, 3=FDD)? ";
MODE = input(prompt);
if MODE == 1
    OMAfreq = SSIcovmu_freq;
    OMAphi = SSIcovmu_modes';
elseif MODE == 2
    OMAfreq = SSIdatFreq;
    OMAphi = SSIdatphi;
elseif MODE == 3
    OMAfreq = FDDmu_freq;
    OMAphi = FDDmu_modes';
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

S = [{'height','numphi1','numphi2','numphi3','numphi4','numphi5','OMAphi1','OMAphi2','OMAphi3','OMAphi4','OMAphi5','OMAfreq'};num2cell(x'),num2cell(phi),num2cell([zeros(1,length(OMAphi))  ;OMAphi]),num2cell([0;OMAfreq])];
T = array2table([num2cell(x'),num2cell(phi),num2cell([zeros(1,length(OMAphi))  ;OMAphi]),num2cell([0;OMAfreq])]);
T.Properties.VariableNames(1:12) = {'height','numphi1','numphi2','numphi3','numphi4','numphi5','OMAphi1','OMAphi2','OMAphi3','OMAphi4','OMAphi5','OMAfreq'};

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
    writetable(T,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap5_SSIcov_vs_simulated.xlsm')
elseif MODE==2
    disp('Modal Assurance Criterion between Numerical modeshapes and SSIdat  : ')
    disp(strcat(num2str(mac)));
    writecell(S,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap5_SSIdat_vs_simulated.xlsm')
elseif MODE == 3
    disp('Modal Assurance Criterion between Numerical modeshapes and FDD  :' )
    disp(strcat(num2str(mac)));
    writecell(S,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap5_FDD_vs_simulated.xlsm')
end
disp('----------------------------------------------------------------------')
disp(strcat('Mode shape accuracy (MAC),1 : ',num2str(dmac(1)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),2 : ',num2str(dmac(2)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),3 : ',num2str(dmac(3)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),4 : ',num2str(dmac(4)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),5 : ',num2str(dmac(5)*100),'%'));
disp(strcat('Mean mode shape accuracy (MAC): ',num2str(mean(dmac)*100),'%'));
disp('----------------------------------------------------------------------')





