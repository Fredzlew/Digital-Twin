% parameters
clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'),genpath('python'))
% Loading modal parameters from OMA 
% SSI
SSIstab2 = readNPY('stab.npy');
% FDD
FDDPSD = readNPY('FDDPSD.npy');
% definition:
% freq = SSIstab(:,1);
% modelorder = SSIstab(:,2);
% label = SSIstab(:,3);
% damp = SSIstab(:,4);

% Only use frequency up to 20 
k=1;

for i = 1:length(SSIstab2)
 if SSIstab2(i,1) <= 20 && SSIstab2(i,4) <= 0.05 && SSIstab2(i,4) >= 0
     SSIstab(k,:) = SSIstab2(i,:);
     k=k+1;
 end
end

% sortering out from label
m = 1;
n = 1;
s = 1;
q = 1;
b = 1;
for i = 1:length(SSIstab)
    if SSIstab(i,3) == 0 % label = 0
        freq0(m) = SSIstab(i,1);
        modelorder0(m) = SSIstab(i,2);
        damp0(m) = SSIstab(i,4);
        m = m+1;
    elseif SSIstab(i,3) == 1 % label = 1
        freq1(n) = SSIstab(i,1);
        modelorder1(n) = SSIstab(i,2);
        damp1(n) = SSIstab(i,4);
        n = n+1;
    elseif SSIstab(i,3) == 2 % label = 2
        freq2(s) = SSIstab(i,1);
        modelorder2(s) = SSIstab(i,2);
        damp2(s) = SSIstab(i,4);
        s = s+1;
    elseif SSIstab(i,3) == 3 % label = 3
        freq3(q) = SSIstab(i,1);
        modelorder3(q) = SSIstab(i,2);
        damp3(q) = SSIstab(i,4);
        q = q+1;    
    elseif SSIstab(i,3) == 4 % label = 4
        freq4(b) = SSIstab(i,1);
        modelorder4(b) = SSIstab(i,2);
        damp4(b) = SSIstab(i,4);
        b = b+1;
    end
end

% plotting stabilization diagram in 3D
figure
plot3(freq0,damp0,modelorder0,'ro',freq1,damp1,modelorder1,'mo',freq2,damp2,modelorder2,'yo',freq3,damp3,modelorder3,'bo',freq4,damp4,modelorder4,'go')
legend('Unstable pole','Stable for frequency','Stable for frequency and damping','Stable for frequency and mode shape','Stable pole','FontSize', 12)
xlabel('Frequency [Hz]','FontSize', 14)
ylabel('Damping ratio [-]','FontSize', 14)
zlabel('Model order','FontSize', 14)
ylim([0 0.05])

%% Plotting stabilization diagram 2D with FDD
close all
% Taking the PSD out of the 3D matrix
for i = 1:length(FDDPSD)
    PSD1 = diag(FDDPSD(:,:,i));
    PSD(i,1) = PSD1(1);
    PSD(i,2) = PSD1(2);
    PSD(i,3) = PSD1(3);
    PSD(i,4) = PSD1(4);
    PSD(i,5) = PSD1(5);
end
% load frequencies
FDDPSDfreq = readNPY('FDDPSDfreq.npy');
f_ = FDDPSDfreq(1:length(FDDPSD));
%plotting all figure together
figure 
hold on
% plot PSD
for i = 1:5
plot(f_,10*log10(PSD(:,i))+56)
end
% plot SSI
plot(freq0,modelorder0,'ro')
plot(freq1,modelorder1,'mo')
plot(freq2,modelorder2,'yo')
plot(freq3,modelorder3,'bo')
plot(freq4,modelorder4,'go')
legend('Unstable pole','Stable for frequency','Stable for frequency and damping','Stable for frequency and mode shape','Stable pole','FontSize', 12)
xlabel('Frequency [Hz]','FontSize', 14)
ylabel('Model order','FontSize', 14)
xlim([0 20])
ylim([0 130])