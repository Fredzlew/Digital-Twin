clc; clear; close all;
%Model Parameters and excitation
%--------------------------------------------------------------------------
% Set global random seed
rng(2)

% Known system matrices
M=[1 0; 0 1];
K=[2 -1; -1 1]*5;
C=0.1*M+0.1*K;
f=2*randn(2,10000);
fs=100;

%Apply modal superposition to get response
%--------------------------------------------------------------------------

n=size(f,1);
dt=1/fs; %sampling rate
[Vectors, Values]=eig(K,M);
Freq=sqrt(diag(Values))/(2*pi); % undamped natural frequency
steps=size(f,2);

Mn=diag(Vectors'*M*Vectors); % uncoupled mass
Cn=diag(Vectors'*C*Vectors); % uncoupled damping
Kn=diag(Vectors'*K*Vectors); % uncoupled stifness
wn=sqrt(diag(Values));
zeta=Cn./(sqrt(2.*Mn.*Kn));  % damping ratio
wd=wn.*sqrt(1-zeta.^2);

fn=Vectors'*f; % generalized input force matrix

t=[0:dt:dt*steps-dt];

for i=1:1:n
    
    h(i,:)=(1/(Mn(i)*wd(i))).*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t); %transfer function of displacement
    hd(i,:)=(1/(Mn(i)*wd(i))).*(-zeta(i).*wn(i).*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)+wd(i).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t)); %transfer function of velocity
    hdd(i,:)=(1/(Mn(i)*wd(i))).*((zeta(i).*wn(i))^2.*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)-zeta(i).*wn(i).*wd(i).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t)-wd(i).*((zeta(i).*wn(i)).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t))-wd(i)^2.*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)); %transfer function of acceleration
    
    qq=conv(fn(i,:),h(i,:))*dt;
    qqd=conv(fn(i,:),hd(i,:))*dt;
    qqdd=conv(fn(i,:),hdd(i,:))*dt;
    
    q(i,:)=qq(1:steps); % modal displacement
    qd(i,:)=qqd(1:steps); % modal velocity
    qdd(i,:)=qqdd(1:steps); % modal acceleration
       
end

x=Vectors*q; %displacement
v=Vectors*qd; %vecloity
a=Vectors*qdd; %vecloity

%Add noise to excitation and response
%--------------------------------------------------------------------------
f2=f+0.05*randn(2,10000);
a2=a+0.05*randn(2,10000);
v2=v+0.05*randn(2,10000);
x2=x+0.05*randn(2,10000);

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
output=x2;
ncols=5000;
nrows=600;
cut=4;  %Identify 2 modes
[Result]=SSID(output,fs,ncols,nrows,cut);    %SSI

%Plot real and identified first modes to compare between them
%--------------------------------------------------------------------------
figure;
plot([0 ; -Vectors(:,1)],[0 1 2],'r*-');
hold on
plot([0  ;Result.Parameters.ModeShape(:,1)],[0 1 2],'go-.');
hold on
plot([0 ; -Vectors(:,2)],[0 1 2],'b^--');
hold on
plot([0  ;Result.Parameters.ModeShape(:,2)],[0 1 2],'mv-');
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