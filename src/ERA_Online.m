clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
%Model Parameters and excitation
%--------------------------------------------------------------------------
% Choose stiffness matrix
prompttt = "Ricky and Johan's or Jan's stiffness matrix: (1=Ricky&Johan, 2=Jan)? ";
prop = input(prompttt);
if prop == 1
    filename = load('modelprop.mat'); % Loads mass and stiffness matrices
elseif prop == 2
    filename = load('modelprop_jan.mat'); % Loads mass and stiffness matrices
end

% Choose of data
prompt = "Use ERA for measured or simulated data (1=measured, 2=simulated(IRF), 3=simulated(Newmark))? ";
ERAdata = input(prompt);
if ERAdata == 1
    % Measurements
    data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
    fss = data(2:6,1:10000)/1000; % Converting mm to m
    f = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]; % Swap columns due to sensor
elseif ERAdata == 2 && prop == 1
    % Simulated data
    data_sim = load('data_sim.mat');
    f = data_sim.dis(:,1:10000);
elseif ERAdata == 2 && prop == 2
    % Simulated data
    data_sim = load('data_sim_jan.mat');
    f = data_sim.dis(:,1:10000);
elseif ERAdata == 3 && prop == 1
    % Simulated data
    data_sim = load('data_sim_newmark.mat');
    f = data_sim.dis_new(:,1:10000);
elseif ERAdata == 3 && prop == 2
    % Simulated data
    data_sim = load('data_sim_newmark_jan.mat');
    f = data_sim.dis_new(:,1:10000);
end


omega_min = 1.70; % Minimum natural frequency
zeta_min = 0.015; % Minimum threshold for desired damping ratio
alpha = zeta_min*omega_min; % Rayleigh damping coefficient
beta = zeta_min/omega_min; % Rayleigh damping coefficient
M=filename.M; % Mass matrix
K=filename.K; % Stiffness matrix
C=0;%alpha*M+beta*K; % Damping matrix using Rayleigh damping
fs=100; % Sampling frequency (1/dt)
n=size(f,1);
dt=1/fs; %sampling rate
% Solve eigenvalue problem to find numerical modal parameters
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

%Identify modal parameters using displacement with added uncertainty
%--------------------------------------------------------------------------
nm = 5; %Number of modes
Y=f; %Displacements
nrows=5*512;%112;%50*(2*nm/5)+1; Best at 2000, realistic around 600
ncols=7400;%7044;%4/5*size(f,2)-nrows-3;    %more than 2/3 of No. of data
inputs=1;     
cut=2*nm;        %Identify 5 modes
shift=10;      %Adjust EMAC sensitivity
EMAC_option=1; %EMAC is calculated only from observability matrix

[Result]=ERA(Y,fs,ncols,nrows,inputs,cut,shift,EMAC_option);  %ERA

% Normalizing OMA mode shapes
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


% plotting the mode shapes
x = [0 filename.H];
phi = [zeros(1,length(Vectors)); Vectors];
fig = figure;
fig.Position=[100 100 1600 700];
for i=1:nm
    subplot(1,nm,i)
    hold on
    plot(phi(:,i),x,'-m')
    if phi(2,i)*phi_ERA(1,i) < 0 % Swap sign on mode shape
        plot([0  ;-phi_ERA(:,i)],x,'go-.');
        plot(-phi_ERA(1:end,i),x(2:end),'g.','markersize',30)
    else
        plot([0  ;phi_ERA(:,i)],x,'go-.');
        plot(phi_ERA(1:end,i),x(2:end),'g.','markersize',30)
    end
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

if ERAdata == 1
    save('.\data\ERAmodal.mat','phi_ERA','ERAFreq');
elseif ERAdata == 2
    save('.\data\ERAmodalsim.mat','phi_ERA','ERAFreq');
elseif ERAdata == 3
    save('.\data\ERAmodalsim_newmark.mat','phi_ERA','ERAFreq');
end