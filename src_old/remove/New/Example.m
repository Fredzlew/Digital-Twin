clc; clear; close all;
%Model Parameters and excitation
%--------------------------------------------------------------------------

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
f2=f+0.01*randn(2,10000);
a2=a+0.01*randn(2,10000);
v2=v+0.01*randn(2,10000);
x2=x+0.01*randn(2,10000);

%Identify modal parameters using displacement with added uncertainty
%--------------------------------------------------------------------------
input=f2;
output=x2;

cut=4;  %Identify 2 modes
[Result1]=SSI(output,fs,cut);         %SSI
[Result2]=DSI(output,input,fs,cut);   %DSI
[Result3]=DSSI(output,input,fs,cut);  %DSSI

%Plot real and identified first modes to compare between them
%--------------------------------------------------------------------------
figure;
plot([0 ; -Vectors(:,1)],[0 1 2],'b*-');
hold on
plot([0  ;Result1.Parameters.ModeShape(:,1)],[0 1 2],'k*-');
hold on
plot([0 ; Result2.Parameters.ModeShape(:,1)],[0 1 2],'g*-');
hold on
plot([0 ; Result3.Parameters.ModeShape(:,1)],[0 1 2],'r*-');
hold off
title('Real and Identified First Mode Shape');
legend('Real','Identified using SSI','Identified using DSI','Identified using DSSI');
xlabel('Amplitude');
ylabel('Floor');
grid on;

%Display real and Identified natural frequencies and damping ratios
%--------------------------------------------------------------------------
disp('Real and Identified Natural Drequencies and Damping Ratios of the First Mode'); disp('');
disp(strcat('Real: Frequency=',num2str(Freq(1)),'Hz',' Damping Ratio=',num2str(zeta(1)*100),'%'));
disp(strcat('SSI: Frequency=',num2str(Result1.Parameters.NaFreq(1)),'Hz',' Damping Ratio=',num2str(Result1.Parameters.DampRatio(1)),'%'));
disp(strcat('DSI: Frequency=',num2str(Result2.Parameters.NaFreq(1)),'Hz',' Damping Ratio=',num2str(Result2.Parameters.DampRatio(1)),'%'));
disp(strcat('DSSI: Frequency=',num2str(Result3.Parameters.NaFreq(1)),'Hz',' Damping Ratio=',num2str(Result3.Parameters.DampRatio(1)),'%'));
