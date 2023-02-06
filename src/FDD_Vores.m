clc;
clear all;
close all;

dt = 0.01; % here you specify the sampling interval
Nf1 = 2^12; % here you specify the number of frequency lines (a radix 2 number, i.e., Nf1=2^n, where n is a natural number)

% below you specify the 1st three natural frequencies identified from the (frequency domain) plots
Fn1e = 1.65;   % 1st natural frequency
Fn2e = 5.03;   % 2nd natural frequency
Fn3e = 7.90;   % 3rd natural frequency
Fn4e = 10.12;  % 4th natural frequency
Fn5e = 11.60;  % 5th natural frequency

%######## from this point onwards do not edit the script ##################
% addpath('.\OMAtoolsA1_5\Signal processing');
data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
ys = (data(2:6,:)/1000)'; % Converting mm to m
y = [ys(:,5),ys(:,4),ys(:,3),ys(:,2),ys(:,1)]; % Swap columns due to sensor

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
    ModeShapes(:,i) = Uf(:,1);
    Sv(:,i) = diag(Sf);
end

% plot singular values
figure(1);
plot(f,10*log(Sv(1:5,:).^0.5));
xlabel('Frequency [Hz]');
ylabel('Sigular Values [dB to unit]');
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
% Loading modelprop for real mode shapes
filename = load('modelprop.mat'); % Loads mass and stiffness matrices

% plotting the mode shapes
x = [0 filename.H];
phi = [zeros(1,length(filename.U)); filename.U];
MSV = [-Vn1,-Vn2,-Vn3,-Vn4,-Vn5];
fig = figure;
fig.Position=[100 100 1600 700];
for i=1:length(Uf)
    subplot(1,length(Uf),i)
    hold on
    plot(phi(:,i),x,'-m')
    plot([0  ;MSV(:,i)],x,'go-.');
    plot(phi(2:end,i),x(2:end),'b.','markersize',30)
%     title(['f = ' num2str(fn(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
    xline(0.0,'--')
    xlim([-1.1,1.1])
    ylim([0,x(end)])
    if i==1
        legend('Numerical','Approximation','Location','northwest')
    else
    end
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);
