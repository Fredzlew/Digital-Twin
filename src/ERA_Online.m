clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
%Model Parameters and excitation
%--------------------------------------------------------------------------
%rng(1) % Set global random seed
data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
fss = data(2:6,1:10000)/1000; % Converting mm to m
f = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]; % Swap columns due to sensor
filename = load('modelprop.mat'); % Loads mass and stiffness matrices
omega_min = 1.70; % Minimum natural frequency
zeta_min = 0.015; % Minimum threshold for desired damping ratio
alpha = zeta_min*omega_min; % Rayleigh damping coefficient
beta = zeta_min/omega_min; % Rayleigh damping coefficient

M=filename.M; % Mass matrix
K=filename.K; % Stiffness matrix
C=0;%alpha*M+beta*K; % Damping matrix using Rayleigh damping
%f=2*randn(5,10000); %Generating 10000 measurements from 2 floors
fs=100; % Sampling frequency (1/dt)

%Apply modal superposition to get response

n=size(f,1);
dt=1/fs; %sampling rate
[Us, Values]=eig(K,M);
Freq=sqrt(diag(Values))/(2*pi); % undamped natural frequency
steps=size(f,2);

% normalizing mode shapes
MVec_x = max(Us); % start normalization
mVec_x = min(Us);
for j = 1:5
    if abs(MVec_x(j)) > abs(mVec_x(j))
        mxVec_x(j) = MVec_x(j);
    else
        mxVec_x(j) = mVec_x(j);
    end
    for l = 1:5
        Vectors(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization

% Mn=diag(Vectors'*M*Vectors); % uncoupled mass
% Cn=diag(Vectors'*C*Vectors); % uncoupled damping
% Kn=diag(Vectors'*K*Vectors); % uncoupled stifness
% wn=sqrt(diag(Values));
% zeta=Cn./(sqrt(2.*Mn.*Kn));  % damping ratio
% wd=wn.*sqrt(1-zeta.^2);
% 
% fn=Vectors'*f; % generalized input force matrix
% 
% t=[0:dt:dt*steps-dt];
% 
% for i=1:1:n
%     
%     h(i,:)=(1/(Mn(i)*wd(i))).*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t); %transfer function of displacement
%     hd(i,:)=(1/(Mn(i)*wd(i))).*(-zeta(i).*wn(i).*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)+wd(i).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t)); %transfer function of velocity
%     hdd(i,:)=(1/(Mn(i)*wd(i))).*((zeta(i).*wn(i))^2.*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)-zeta(i).*wn(i).*wd(i).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t)-wd(i).*((zeta(i).*wn(i)).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t))-wd(i)^2.*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)); %transfer function of acceleration
%     
%     qq=conv(fn(i,:),h(i,:))*dt;
%     qqd=conv(fn(i,:),hd(i,:))*dt;
%     qqdd=conv(fn(i,:),hdd(i,:))*dt;
%     
%     q(i,:)=qq(1:steps); % modal displacement
%     qd(i,:)=qqd(1:steps); % modal velocity
%     qdd(i,:)=qqdd(1:steps); % modal acceleration
%        
% end
% 
% x=Vectors*q; %displacement
% v=Vectors*qd; %vecloity
% a=Vectors*qdd; %vecloity
% 
% %Add noise to excitation and response
% %--------------------------------------------------------------------------
% f2=f;%+0.01*randn(2,10000);
% a2=a;%+0.01*randn(2,10000);
% v2=v;%+0.01*randn(2,10000);
% x2=x;%+0.01*randn(2,10000);

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
nm = 5; %Number of modes
Y=f; %Displacements
ncols=4/5*length(f);    %more than 2/3 of No. of data
nrows=20*2*nm/5+1;     %more than 20 * number of modes
inputs=1;     
cut=2*nm;        %Identify 5 modes
shift=10;      %Adjust EMAC sensitivity
EMAC_option=1; %EMAC is calculated only from observability matrix

[Result]=ERA(Y,fs,ncols,nrows,inputs,cut,shift,EMAC_option);  %ERA

% normalizing SSI mode shapes
Us = Result.Parameters.ModeShape;
MVec_x = max(Us); % start normalization
mVec_x = min(Us);
for j = 1:nm
    if abs(MVec_x(j)) > abs(mVec_x(j))
        mxVec_x(j) = MVec_x(j);
    else
        mxVec_x(j) = mVec_x(j);
    end
    for l = 1:5
        phi_ERA(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization

%Plot real and identified first modes to compare between them
%--------------------------------------------------------------------------
% figure;
% plot([0 ; Vectors(:,1)],[0 1 2 3 4 5],'b*-');
% hold on
% plot([0  ;Result.Parameters.ModeShape(:,1)],[0 1 2 3 4 5],'r*-.');
% hold on
% plot([0 ; -Vectors(:,2)],[0 1 2 3 4 5],'g*--');
% hold on
% plot([0  ;Result.Parameters.ModeShape(:,3)],[0 1 2 3 4 5],'k*:');
% hold off
% title('Real and Identified Mode Shapes');
% legend('Mode 1 (Real)','Mode 1 (Identified using ERA)','Mode 2 (Real)','Mode 2 (Identified using ERA)');
% xlabel('Amplitude');
% ylabel('Floor');
% grid on;

% plotting the mode shapes
x = [0 filename.H];
phi = [zeros(1,length(Vectors)); Vectors];
fig = figure;
fig.Position=[100 100 1600 700];
for i=1:nm
    subplot(1,nm,i)
    hold on
    if i == 4
        plot(phi(:,i),x,'-m')
        plot([0  ;-phi_ERA(:,i)],x,'go-.');
        plot(phi(2:end,i),x(2:end),'b.','markersize',30)
        title(['f = ' num2str(Result.Parameters.NaFreq(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
        xline(0.0,'--')
        xlim([-1.1,1.1])
        ylim([0,x(end)])
    else
        plot(phi(:,i),x,'-m')
        plot([0  ;phi_ERA(:,i)],x,'go-.');
        plot(phi(2:end,i),x(2:end),'b.','markersize',30)
        title(['f = ' num2str(Result.Parameters.NaFreq(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
        xline(0.0,'--')
        xlim([-1.1,1.1])
        ylim([0,x(end)])
        if i==1
            legend('Numerical','Approximation','Location','northwest')
        else
        end
    end
end
sgtitle('Numerical vs ERA Online','FontSize',20) 
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);
ERAFreq = Result.Parameters.NaFreq;
%Display real and Identified natural frequencies and damping ratios
%--------------------------------------------------------------------------
disp('Real and Identified Natural Drequencies and Damping Ratios of the First Mode'); disp('');
%disp(strcat('Real: Frequency=',num2str(Freq(1)),'Hz',' Damping Ratio=',num2str(zeta(1)*100),'%'));
disp(strcat('ERA: Frequency=',num2str(Result.Parameters.NaFreq(1)),'Hz',' Damping Ratio=',num2str(Result.Parameters.DampRatio(1)),'%'));
disp(strcat('CMI=',num2str(Result.Indicators.CMI(1))));

disp('Real and Identified Natural Drequencies and Damping Ratios of the Second Mode'); disp('');
%disp(strcat('Real: Frequency=',num2str(Freq(2)),'Hz',' Damping Ratio=',num2str(zeta(2)*100),'%'));
disp(strcat('ERA: Frequency=',num2str(Result.Parameters.NaFreq(2)),'Hz',' Damping Ratio=',num2str(Result.Parameters.DampRatio(2)),'%'));
disp(strcat('CMI=',num2str(Result.Indicators.CMI(2))));

save('.\data\ERAmodal.mat','phi_ERA','ERAFreq');