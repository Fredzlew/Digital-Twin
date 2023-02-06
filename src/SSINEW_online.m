clc; clear; close all;
%Model Parameters and excitation
%--------------------------------------------------------------------------

data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
f = data(2:6,1:10000)/1000; % Converting mm to m
%f = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]; % Swap columns due to sensor
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
f2=f; %+0.01*randn(2,10000);
a2=a;%+0.01*randn(2,10000);
v2=v;%+0.01*randn(2,10000);
x2=x;%+0.01*randn(2,10000);

%Identify modal parameters using displacement with added uncertainty
%--------------------------------------------------------------------------
input=f2;
output=x2;
nm = 5; % number of modes

cut=nm*2;  %Identify 2 modes
[Result]=SSI(output,fs,cut);         %SSI
[Result2]=DSI(output,input,fs,cut);   %DSI
[Result3]=DSSI(output,input,fs,cut);  %DSSI

% normalizing SSI mode shapes
Us = Result.Parameters.ModeShape;
MVec_x = max(Us); % start normalization
mVec_x = min(Us);
for j = 1:5
    if abs(MVec_x(j)) > abs(mVec_x(j))
        mxVec_x(j) = MVec_x(j);
    else
        mxVec_x(j) = mVec_x(j);
    end
    for l = 1:5
        phi_SSI(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization


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
    if i == 5
        plot(phi(:,i),x,'-m')
        plot([0  ;-phi_SSI(:,i)],x,'go-.');
        plot(phi(2:end,i),x(2:end),'b.','markersize',30)
        title(['f = ' num2str(Result.Parameters.NaFreq(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
        xline(0.0,'--')
        xlim([-1.1,1.1])
        ylim([0,x(end)])
    else
        plot(phi(:,i),x,'-m')
        plot([0  ;phi_SSI(:,i)],x,'go-.');
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
sgtitle('Numerical vs SSI Online','FontSize',20) 
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);

% figure;
% plot([0 ; Vectors(:,1)],[0 1 2 3 4 5],'b*-');
% hold on
% plot([0  ;Result1.Parameters.ModeShape(:,1)],[0 1 2 3 4 5],'k*-');
% hold on
% plot([0 ; Result2.Parameters.ModeShape(:,1)],[0 1 2 3 4 5],'g*-');
% hold on
% plot([0 ; Result3.Parameters.ModeShape(:,1)],[0 1 2 3 4 5],'r*-');
% hold off
% title('Real and Identified First Mode Shape');
% legend('Real','Identified using SSI','Identified using DSI','Identified using DSSI');
% xlabel('Amplitude');
% ylabel('Floor');
% grid on;

%Display real and Identified natural frequencies and damping ratios
%--------------------------------------------------------------------------
disp('Real and Identified Natural Drequencies and Damping Ratios of the First Mode'); disp('');
disp(strcat('Real: Frequency=',num2str(Freq(1)),'Hz',' Damping Ratio=',num2str(zeta(1)*100),'%'));
disp(strcat('SSI: Frequency=',num2str(Result.Parameters.NaFreq(1)),'Hz',' Damping Ratio=',num2str(Result.Parameters.DampRatio(1)),'%'));
disp(strcat('DSI: Frequency=',num2str(Result2.Parameters.NaFreq(1)),'Hz',' Damping Ratio=',num2str(Result2.Parameters.DampRatio(1)),'%'));
disp(strcat('DSSI: Frequency=',num2str(Result3.Parameters.NaFreq(1)),'Hz',' Damping Ratio=',num2str(Result3.Parameters.DampRatio(1)),'%'));
