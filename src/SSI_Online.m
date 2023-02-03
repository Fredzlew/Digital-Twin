clc; clear; close all;
%Model Parameters and excitation
%--------------------------------------------------------------------------
rng(1) % Set global random seed
addpath('C:\Users\User\Danmarks Tekniske Universitet\Frederik Emil Serritzlew - Kandidat\Experimental_data\Anela')
FullInFName = 'data_1_2_1.txt'; % Loading displacement data
f = dlmread(FullInFName,'',7,1)'; % Read data from file
f = f(1:5,1:10000)/1000; % Converting mm to m
filename = load('modelprop.mat'); % Loads mass and stiffness matrices

M=filename.M; % Mass matrix
K=filename.K; % Stiffness matrix
C=0.1*M+0.1*(K/1000); % Damping matrix using Rayleigh damping
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
figure;
subplot(3,2,1)
plot(t,f(1,:)); xlabel('Time (sec)');  ylabel('Force1'); title('First Floor');
subplot(3,2,2)
plot(t,f(2,:)); xlabel('Time (sec)');  ylabel('Force2'); title('Second Floor');
subplot(3,2,3)
plot(t,x(1,:)); xlabel('Time (sec)');  ylabel('DSP1');
subplot(3,2,4)
plot(t,x(2,:)); xlabel('Time (sec)');  ylabel('DSP2');
subplot(3,2,5)
plot(t,x2(1,:)); xlabel('Time (sec)');  ylabel('DSP1+Noise');
subplot(3,2,6)
plot(t,x2(2,:)); xlabel('Time (sec)');  ylabel('DSP2+Noise');

%Identify modal parameters using displacement with added uncertainty
%--------------------------------------------------------------------------
output=x2; % Displacements
ncols=7000; % More than 2/3*number of samples
nrows=600; % More than 20*number of sensors
cut=10;  % cut=4 -> 2 modes, cut=10 -> 5 modes
[Result]=SSID(output,fs,ncols,nrows,cut);    %SSI

%Plot real and identified first modes to compare between them
%--------------------------------------------------------------------------
figure;
plot([0 ; -Vectors(:,1)],'r*-');
hold on
plot([0  ;Result.Parameters.ModeShape(:,1)],'go-.');
hold on
plot([0 ; -Vectors(:,2)],'b^--');
hold on
plot([0  ;Result.Parameters.ModeShape(:,2)],'mv-');
hold off
title('Real and Identified Mode Shapes');
legend('Mode 1 (Real)','Mode 1 (Identified using SSI)','Mode 2 (Real)','Mode 2 (Identified using SSI)');
xlabel('Amplitude');
ylabel('Floor');
grid on;
daspect([1 1 1]);

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