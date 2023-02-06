clc; clear; close all;
%Model Parameters and excitation
%--------------------------------------------------------------------------
%rng(1) % Set global random seed
data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
f = data(2:6,1:10000)/1000; % Converting mm to m
filename = load('modelprop.mat'); % Loads mass and stiffness matrices
omega_min = 1.70; % Minimum natural frequency
zeta_min = 0.015; % Minimum threshold for desired damping ratio
alpha = zeta_min*omega_min; % Rayleigh damping coefficient
beta = zeta_min/omega_min; % Rayleigh damping coefficient

M=filename.M; % Mass matrix
K=filename.K; % Stiffness matrix
C=alpha*M+beta*K; % Damping matrix using Rayleigh damping
%f=2*randn(2,10000); %Generating 10000 measurements from 2 floors
fs=100; % Sampling frequency (1/dt)

%Apply modal superposition to get response
%--------------------------------------------------------------------------

n=size(f,1); % Number of floors/sensors
dt=1/fs; %sampling rate
[Vectors, Values]=eig(K,M); % Solving eigenvalueproblem (undamped)
Freq=sqrt(diag(Values))/(2*pi); % Undamped natural frequency
steps=size(f,2); % Number of samples

% Then the natural frequencies, including damping are calculated
Mn=diag(Vectors'*M*Vectors); % uncoupled mass
Cn=diag(Vectors'*C*Vectors); % uncoupled damping
Kn=diag(Vectors'*K*Vectors); % uncoupled stifness
wn=sqrt(diag(Values)); % Undamped natural frequencies
zeta=Cn./(sqrt(2.*Mn.*Kn));  % damping ratio
wd=wn.*sqrt(1-zeta.^2); % Damped natural frequencies

fn=Vectors'*f; % generalized input force matrix (Modal load)

t=[0:dt:dt*steps-dt]; % Time for each measurement

for i=1:1:n
    
    h(i,:)=(1/(Mn(i)*wd(i))).*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t); %transfer function of displacement
    hd(i,:)=(1/(Mn(i)*wd(i))).*(-zeta(i).*wn(i).*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)+wd(i).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t)); %transfer function of velocity
    hdd(i,:)=(1/(Mn(i)*wd(i))).*((zeta(i).*wn(i))^2.*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)-zeta(i).*wn(i).*wd(i).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t)-wd(i).*((zeta(i).*wn(i)).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t))-wd(i)^2.*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)); %transfer function of acceleration
    
    qq=conv(fn(i,:),h(i,:))*dt; % modal displacement
    qqd=conv(fn(i,:),hd(i,:))*dt; % modal velocity
    qqdd=conv(fn(i,:),hdd(i,:))*dt; % modal acceleration
    
    q(i,:)=qq(1:steps); % Picking the displacements within the time frame
    qd(i,:)=qqd(1:steps); % Picking the velocities within the time frame
    qdd(i,:)=qqdd(1:steps); % Picking the accelerations within the time frame
       
end

x=Vectors*q; % Real displacements
v=Vectors*qd; % Real velocities
a=Vectors*qdd; % Real accelerations

%Add noise to excitation and response
%--------------------------------------------------------------------------
f2=f;%+0.05*randn(2,10000);
a2=a;%+0.05*randn(2,10000);
v2=v;%+0.05*randn(2,10000);
x2=x;%+0.05*randn(2,10000);

%Plot displacement of first floor without and with noise
%--------------------------------------------------------------------------
% figure;
% subplot(3,2,1)
% plot(t,f(1,:)); xlabel('Time (sec)');  ylabel('Force1'); title('First Floor');
% subplot(3,2,2)
% plot(t,f(2,:)); xlabel('Time (sec)');  ylabel('Force2'); title('Second Floor');
% subplot(3,2,3)
% plot(t,x(1,:)); xlabel('Time (sec)');  ylabel('DSP1');
% subplot(3,2,4)
% plot(t,x(2,:)); xlabel('Time (sec)');  ylabel('DSP2');
% subplot(3,2,5)
% plot(t,x2(1,:)); xlabel('Time (sec)');  ylabel('DSP1+Noise');
% subplot(3,2,6)
% plot(t,x2(2,:)); xlabel('Time (sec)');  ylabel('DSP2+Noise');

%Identify modal parameters using displacement with added uncertainty
%--------------------------------------------------------------------------
nm = 5; %number of modes
output=x2; % Displacements
ncols=4/5*length(f); % More than 2/3*number of samples
nrows=50*nm; % More than 20*number of sensors
cut=2*nm;  % cut=4 -> 2 modes, cut=10 -> 5 modes
[Result]=SSID(output,fs,ncols,nrows,cut);    %SSI

%Plot real and identified first modes to compare between them
%--------------------------------------------------------------------------
% plotting the mode shapes
x = [0 filename.H];
phi = [zeros(1,length(Vectors)); Vectors];
fig = figure;
fig.Position=[100 100 1600 700];
for i=1:nm
    subplot(1,nm,i)
    hold on
    plot(phi(:,i),x,'-m')
    plot([0  ;Result.Parameters.ModeShape(:,i)],x,'go-.');
    plot(phi(2:end,i),x(2:end),'b.','markersize',30)
%     title(['f = ' num2str(fn(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
    xline(0.0,'--')
    xlim([-1.1,1.1])
    ylim([0,x(end)])
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);

% figure;
% plot([0 ; Vectors(:,1)],[0 1 2 3 4 5],'r*-');
% hold on
% plot([0  ;Result.Parameters.ModeShape(:,1)],[0 1 2 3 4 5],'go-.');
% hold on
% plot([0 ; -Vectors(:,2)],[0 1 2 3 4 5],'b^--');
% hold on
% plot([0  ;Result.Parameters.ModeShape(:,2)],[0 1 2 3 4 5],'mv-');
% hold off
% title('Real and Identified Mode Shapes');
% legend('Mode 1 (Real)','Mode 1 (Identified using SSI)','Mode 2 (Real)','Mode 2 (Identified using SSI)');
% xlabel('Amplitude');
% ylabel('Floor');
% grid on;
% daspect([1 1 1]);

%Display real and Identified natural frequencies and damping ratios
%--------------------------------------------------------------------------
disp('Real and Identified Natural Drequencies and Damping Ratios of the First Mode'); 
disp(strcat('Real: Frequency=',num2str(Freq(1)),'Hz',' Damping Ratio=',num2str(zeta(1)*100),'%'));
disp(strcat('SSI: Frequency=',num2str(Result.Parameters.NaFreq(1)),'Hz',' Damping Ratio=',num2str(Result.Parameters.DampRatio(1)),'%'));
disp(strcat('CMI of The Identified Mode=',num2str(Result.Indicators.CMI(1)),'%'));
disp('-----------')
disp('Real and Identified Natural Drequencies and Damping Ratios of the Second Mode');
disp(strcat('Real: Frequency=',num2str(Freq(2)),'Hz',' Damping Ratio=',num2str(zeta(2)*100),'%'));
disp(strcat('SSI: Frequency=',num2str(Result.Parameters.NaFreq(2)),'Hz',' Damping Ratio=',num2str(Result.Parameters.DampRatio(2)),'%'));
disp(strcat('CMI of The Identified Mode=',num2str(Result.Indicators.CMI(2)),'%'));