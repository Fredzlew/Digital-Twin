clc;
clear all;
close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
Nf1 = 2^14; % here you specify the number of frequency lines (a radix 2 number, i.e., Nf1=2^n, where n is a natural number)


% data from numerical model
prompttt = "Ricky and Johan's or Jan's stiffness matrix: (1=Ricky&Johan, 2=Jan)? ";
prop = input(prompttt);
if prop == 1
    filename = load('modelprop.mat'); % Loads mass and stiffness matrices
elseif prop == 2
    filename = load('modelprop_jan.mat'); % Loads mass and stiffness matrices
end
% the numerical model frequencies
fnumerical = filename.fn;

%######## from this point onwards do not edit the script ##################
% addpath('.\OMAtoolsA1_5\Signal processing');
% Measured data and Simulated data
prompt = "Use FDD for measured or simulated data (1=measured, 2=simulated(IRF), 3=simulated(Newmark))? ";
FDDdata = input(prompt);
if FDDdata == 1
    % below you specify the 1st three natural frequencies identified from the (frequency domain) plots
    % Frequencies with measured data
    Fn1e = 1.65;   % 1st natural frequency
    Fn2e = 5.03;   % 2nd natural frequency
    Fn3e = 7.90;   % 3rd natural frequency
    Fn4e = 10.12;  % 4th natural frequency
    Fn5e = 11.60;  % 5th natural frequency
    % measured data
    data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
    ys = (data(2:6,:)/1000)'; % Converting mm to m
    y = [ys(:,5),ys(:,4),ys(:,3),ys(:,2),ys(:,1)]; % Swap columns due to sensor
    dt = 0.01; % here you specify the sampling interval
elseif FDDdata == 2 && prop == 1
    % below you specify the 1st three natural frequencies identified from the (frequency domain) plots
    % Frequencies with simulated data
    Fn1e = 1.709;   % 1st natural frequency
    Fn2e = 5.072;   % 2nd natural frequency
    Fn3e = 7.947;   % 3rd natural frequency
    Fn4e = 10.077;  % 4th natural frequency
    Fn5e = 11.408;  % 5th natural frequency
    % Simulated data
    data_sim = load('data_sim.mat');
    y=data_sim.dis';
    dt = 0.01; % here you specify the sampling interval
elseif FDDdata == 2 && prop == 2
    % below you specify the 1st three natural frequencies identified from the (frequency domain) plots
    % Frequencies with simulated data
    Fn1e = 1.746;   % 1st natural frequency
    Fn2e = 5.188;   % 2nd natural frequency
    Fn3e = 8.124;   % 3rd natural frequency
    Fn4e = 10.303;  % 4th natural frequency
    Fn5e = 11.652;  % 5th natural frequency
    % Simulated data
    data_sim = load('data_sim_jan.mat');
    y=data_sim.dis';
    dt = 0.01; % here you specify the sampling interval
elseif FDDdata == 3 && prop == 1
    % below you specify the 1st three natural frequencies identified from the (frequency domain) plots
    % Frequencies with simulated data
    Fn1e = 1.709;   % 1st natural frequency
    Fn2e = 5.066;   % 2nd natural frequency
    Fn3e = 7.935;   % 3rd natural frequency
    Fn4e = 10.071;   % 4th natural frequency
    Fn5e = 11.414;  % 5th natural frequency
    % Simulated data
    data_sim = load('data_sim_newmark.mat');
    y=data_sim.dis_new';
    dt = 0.001; % here you specify the sampling interval
elseif FDDdata == 3 && prop == 2
    % below you specify the 1st three natural frequencies identified from the (frequency domain) plots
    % Frequencies with simulated data
    Fn1e = 1.755;   % 1st natural frequency
    Fn2e = 5.188;   % 2nd natural frequency
    Fn3e = 8.118;   % 3rd natural frequency
    Fn4e = 10.300;  % 4th natural frequency
    Fn5e = 11.643;  % 5th natural frequency
    % Simulated data
    data_sim = load('data_sim_newmark_jan.mat');
    y=data_sim.dis_new';
    dt = 0.001; % here you specify the sampling interval
end

Nf = Nf1 + 1;
Ns = 2 * Nf1;

fs = 1/dt; % sampling frenquecy
r = 2;     % filter factor (determines frequency range of interest)
fsc = fs/(2*r); % cut-off sampling frequency

[N, No] = size(y); % compute the number of data samples and the number of 
                   % sensors in the vibration data file

yd = fftfilt1(y.', Ns, dt, 'detrend', 0); % remove trends from accelaration 
                                          % measurements

yf=fftfilt1(yd, Ns, dt, 'decimate', r); % decimate the accelaration 
                                        % measurements
dtr = dt*r; % sampling interval after decimation

[f, Gyy] = fftspec(yf, Ns, dtr); % compute the spectrum matrix

for i=1:Nf
    [Uf,Sf,Vf] = svd(Gyy(:,:,i));
    if i==1
        %disp(Uf)
    end
    ModeShapes(:,i) = Uf(:,1);
    Sv(:,i) = diag(Sf);
end

% plot singular values
figure(1);
plot(f,10*log10(Sv(1:5,:).^0.5));
xlabel('Frequency [Hz]');
ylabel('Singular Values [dB to unit]');
set(gca, 'xgrid', 'on', 'ygrid', 'on');

[~, idF1] = min(abs(f-Fn1e));
[~, idF2] = min(abs(f-Fn2e));
[~, idF3] = min(abs(f-Fn3e));
[~, idF4] = min(abs(f-Fn4e));
[~, idF5] = min(abs(f-Fn5e));

Fn1 = f(idF1);
Fn2 = f(idF2);
Fn3 = f(idF3);
Fn4 = f(idF4);
Fn5 = f(idF5);
fn = [Fn1,Fn2,Fn3,Fn4,Fn5];

% Natural frequnecies

% sprintf('1st Natural Frequency: %0.4f', Fn1)
% sprintf('2nd Natural Frequency: %0.4f', Fn2)
% sprintf('3rd Natural Frequency: %0.4f', Fn3)
% sprintf('4th Natural Frequency: %0.4f', Fn4)
% sprintf('5th Natural Frequency: %0.4f', Fn5)

% Mode shape vectors

Vn1 = abs(ModeShapes(:,idF1)).*cos(angle(ModeShapes(:,idF1)));
Vn2 = abs(ModeShapes(:,idF2)).*cos(angle(ModeShapes(:,idF2)));
Vn3 = abs(ModeShapes(:,idF3)).*cos(angle(ModeShapes(:,idF3)));
Vn4 = abs(ModeShapes(:,idF4)).*cos(angle(ModeShapes(:,idF4)));
Vn5 = abs(ModeShapes(:,idF5)).*cos(angle(ModeShapes(:,idF5)));

% normalizing FDD mode shapes
Us = [Vn1,Vn2,Vn3,Vn4,Vn5];
MVec_x = max(Us); % start normalization
mVec_x = min(Us);
for j = 1:5
    if abs(MVec_x(j)) > abs(mVec_x(j))
        mxVec_x(j) = MVec_x(j);
    else
        mxVec_x(j) = mVec_x(j);
    end
    for l = 1:5
        phi_FDD(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization
% plotting the mode shapes
x = [0 filename.H];
phi = [zeros(1,length(filename.U)); filename.U];
fig = figure;
fig.Position=[100 100 1600 700];
for i=1:length(Uf)
    subplot(1,length(Uf),i)
    hold on
    plot(phi(:,i),x,'-m')
    plot([0  ;phi_FDD(:,i)],x,'go-.');
    plot(phi(2:end,i),x(2:end),'b.','markersize',30)
    title(['f = ' num2str(fn(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
    xline(0.0,'--')
    xlim([-1.1,1.1])
    ylim([0,x(end)])
    if i==1
        legend('Numerical','Approximation','Location','northwest')
    else
    end
end

sgtitle('Numerical vs FDD Vores','FontSize',20) 
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);

if FDDdata == 1
    save('.\data\FDDmodal.mat','phi_FDD','fn');
elseif FDDdata == 2
    save('.\data\FDDmodalsim.mat','phi_FDD','fn');
elseif FDDdata == 3
    save('.\data\FDDmodalsim_newmark.mat','phi_FDD','fn');
end

disp(filename.fn)
disp(fn')