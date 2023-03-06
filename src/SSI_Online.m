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
prompt = "Use SSI for measured or simulated data (1=measured, 2=simulated(IRF), 3=simulated(Newmark))? ";
SSIdata = input(prompt);
if SSIdata == 1
    % Measurements
    data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
    fss = data(2:6,1:10000)/1000; % Converting mm to m
    f = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]; % Swap columns due to sensor
elseif SSIdata == 2 && prop == 1
    % Simulated data
    data_sim = load('data_sim.mat');
    f = data_sim.dis(:,1:10000);
elseif SSIdata == 2 && prop == 2
    % Simulated data
    data_sim = load('data_sim_jan.mat');
    f = data_sim.dis(:,1:10000);
elseif SSIdata == 3 && prop == 1
    % Simulated data
    data_sim = load('data_sim_newmark.mat');
    f = data_sim.dis_new(:,1:10000);
elseif SSIdata == 3 && prop == 2
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

%Apply modal superposition to get response
%--------------------------------------------------------------------------

n=size(f,1); % Number of floors/sensors
dt=1/fs; %sampling rate
[Us, Values]=eig(K,M); % Solving eigenvalueproblem (undamped)
Freq=sqrt(diag(Values))/(2*pi); % Undamped natural frequency
steps=size(f,2); % Number of samples

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
nm = 5; %number of modes
output=f; % Displacements
ncols=7400;%4/5*length(f); % More than 2/3*number of samples
nrows=5*512;%;%50*nm; % More than 20*number of sensors % Best at 3970, realsitic around 600
cut=2*nm;  % cut=4 -> 2 modes, cut=10 -> 5 modes
[Result]=SSID(output,fs,ncols,nrows,cut);    %SSI

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
    plot(phi(:,i),x,'-m')
    if phi(2,i)*phi_SSI(1,i) < 0 % Swap sign on mode shape
        plot([0  ;-phi_SSI(:,i)],x,'go-.');
        plot(-phi_SSI(1:end,i),x(2:end),'g.','markersize',30)
    else
        plot([0  ;phi_SSI(:,i)],x,'go-.');
        plot(phi_SSI(1:end,i),x(2:end),'g.','markersize',30)
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
sgtitle('Numerical vs SSI Online','FontSize',20) 
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);


SSIFreq = Result.Parameters.NaFreq;
%Display real and Identified natural frequencies and damping ratios
%--------------------------------------------------------------------------
disp('Real and Identified Natural Drequencies and Damping Ratios of the First Mode'); 
%disp(strcat('Real: Frequency=',num2str(Freq(1)),'Hz',' Damping Ratio=',num2str(zeta(1)*100),'%'));
disp(strcat('SSI: Frequency=',num2str(Result.Parameters.NaFreq(1)),'Hz',' Damping Ratio=',num2str(Result.Parameters.DampRatio(1)),'%'));
disp(strcat('CMI of The Identified Mode=',num2str(Result.Indicators.CMI(1)),'%'));
disp('-----------')
disp('Real and Identified Natural Drequencies and Damping Ratios of the Second Mode');
%disp(strcat('Real: Frequency=',num2str(Freq(2)),'Hz',' Damping Ratio=',num2str(zeta(2)*100),'%'));
disp(strcat('SSI: Frequency=',num2str(Result.Parameters.NaFreq(2)),'Hz',' Damping Ratio=',num2str(Result.Parameters.DampRatio(2)),'%'));
disp(strcat('CMI of The Identified Mode=',num2str(Result.Indicators.CMI(2)),'%'));

if SSIdata == 1
    save('.\data\SSImodal.mat','phi_SSI','SSIFreq');
elseif SSIdata == 2
    save('.\data\SSImodalsim.mat','phi_SSI','SSIFreq');
elseif SSIdata == 3
    save('.\data\SSImodalsim_newmark.mat','phi_SSI','SSIFreq');
end


disp(filename.fn)
disp(Result.Parameters.NaFreq)