% parameters
clc; clear; close all;
addpath(genpath('..\data'),genpath('..\npy-matlab-master'))
% Loading modal parameters from OMA 
promptt = "High damping or no damping? (1 = High and 2 = no damp): ";
xx = input(promptt);
if xx == 1
    % SSI
    SSIstab2 = readNPY('..\data\experimental_data\Modal_par\SSIstab_5_2_1.npy');
    % FDD
    FDDPSD = readNPY('..\data\experimental_data\Modal_par\FDDPSD_5_2_1.npy');
    FDDPSDfreq = readNPY('..\data\experimental_data\Modal_par\FDDPSDfreq_5_2_1.npy');
elseif xx == 2
    % SSI
    SSIstab2 = readNPY('..\data\experimental_data\Modal_par\SSIstab_no_damp.npy');
    % FDD
    FDDPSD = readNPY('..\data\experimental_data\Modal_par\FDDPSD_no_damp.npy');
    FDDPSDfreq = readNPY('..\data\experimental_data\Modal_par\FDDPSDfreq_no_damp.npy');
end
% definition:
% freq = SSIstab(:,1);
% modelorder = SSIstab(:,2);
% label = SSIstab(:,3);
% damp = SSIstab(:,4);

% Only use frequency up to 20 
k=1;

for i = 1:length(SSIstab2)
 if SSIstab2(i,1) <= 12 
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
plot3(freq0,damp0,modelorder0*2,'r.',freq1,damp1,modelorder1*2,'m.',freq2,damp2,modelorder2*2,'y.',freq3,damp3,modelorder3*2,'b.',freq4,damp4,modelorder4*2,'g.','MarkerSize',20)
legend('Unstable pole','Stable for frequency','Stable for frequency and damping','Stable for frequency and mode shape','Stable pole','FontSize', 12)
xlabel('Frequency [Hz]','FontSize', 14)
ylabel('Damping ratio [-]','FontSize', 14)
zlabel('Model order','FontSize', 14)
ylim([0 0.01])
zlim([0 100])



%% Plotting stabilization diagram 2D with FDD
close all
% Taking the PSD out of the 3D matrix
for i = 1:length(FDDPSD)
    PSD1 = diag(FDDPSD(:,:,i));
    PSD(i,1) = 10*log10(PSD1(1))+35;
    PSD(i,2) = 10*log10(PSD1(2))+35;
    PSD(i,3) = 10*log10(PSD1(3))+35;
    PSD(i,4) = 10*log10(PSD1(4))+35;
    PSD(i,5) = 10*log10(PSD1(5))+35;
end

% load frequencies
f_ = FDDPSDfreq(1:length(FDDPSD));

%plotting all figure together
figure 
hold on
% plot SSI when plotting remember: 2 times modelorder
scatter(freq0,modelorder0*2,'r','filled')
scatter(freq1,modelorder1*2,'m','filled')
scatter(freq2,modelorder2*2,'y','filled')
scatter(freq3,modelorder3*2,'b','filled')
scatter(freq4,modelorder4*2,'g','filled')
% plot PSD
for i = 1:5
plot(f_,(PSD(:,i)))
end
hold off

legend('Unstable pole','Stable for frequency','Stable for frequency and damping','Stable for frequency and mode shape','Stable pole','FontSize', 12)
xlabel('Frequency [Hz]','FontSize', 14)
ylabel('Model order','FontSize', 14)
xlim([0 12])
ylim([0 100])

if xx == 1
    T = array2table([num2cell(f_),num2cell(PSD)]);
    T1 = array2table([num2cell(freq0'),num2cell((modelorder0*2)'),num2cell(damp0')]);
    T2 = array2table([num2cell(freq1'),num2cell((modelorder1*2)'),num2cell(damp1')]);
    T3 = array2table([num2cell(freq2'),num2cell((modelorder2*2)'),num2cell(damp2')]);
    T4 = array2table([num2cell(freq3'),num2cell((modelorder3*2)'),num2cell(damp3')]);
    T5 = array2table([num2cell(freq4'),num2cell((modelorder4*2)'),num2cell(damp4')]);
    T.Properties.VariableNames(1:6) = {'f_','PSD1','PSD2','PSD3','PSD4','PSD5'};
    T1.Properties.VariableNames(1:3) = {'freq0','modelorder0','damp0'};
    T2.Properties.VariableNames(1:3) = {'freq1','modelorder1','damp1'};
    T3.Properties.VariableNames(1:3) = {'freq2','modelorder2','damp2'};
    T4.Properties.VariableNames(1:3) = {'freq3','modelorder3','damp3'};
    T5.Properties.VariableNames(1:3) = {'freq4','modelorder4','damp4'};
    writetable(T,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_PSD_highdamp.csv','Delimiter',';')
    writetable(T1,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_stab1_highdamp.txt','Delimiter',' ')
    writetable(T2,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_stab2_highdamp.txt','Delimiter',' ')
    writetable(T3,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_stab3_highdamp.txt','Delimiter',' ')
    writetable(T4,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_stab4_highdamp.txt','Delimiter',' ')
    writetable(T5,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_stab5_highdamp.txt','Delimiter',' ')
elseif xx == 2
    T = array2table([num2cell(f_),num2cell(PSD)]);
    T1 = array2table([num2cell(freq0'),num2cell((modelorder0*2)'),num2cell(damp0')]);
    T2 = array2table([num2cell(freq1'),num2cell((modelorder1*2)'),num2cell(damp1')]);
    T3 = array2table([num2cell(freq2'),num2cell((modelorder2*2)'),num2cell(damp2')]);
    T4 = array2table([num2cell(freq3'),num2cell((modelorder3*2)'),num2cell(damp3')]);
    T5 = array2table([num2cell(freq4'),num2cell((modelorder4*2)'),num2cell(damp4')]);
    T.Properties.VariableNames(1:6) = {'f_','PSD1','PSD2','PSD3','PSD4','PSD5'};
    T1.Properties.VariableNames(1:3) = {'freq0','modelorder0','damp0'};
    T2.Properties.VariableNames(1:3) = {'freq1','modelorder1','damp1'};
    T3.Properties.VariableNames(1:3) = {'freq2','modelorder2','damp2'};
    T4.Properties.VariableNames(1:3) = {'freq3','modelorder3','damp3'};
    T5.Properties.VariableNames(1:3) = {'freq4','modelorder4','damp4'};
    writetable(T,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_PSD_noamp.csv','Delimiter',';')
    writetable(T1,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_stab1_nodamp.txt','Delimiter',' ')
    writetable(T2,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_stab2_nodamp.txt','Delimiter',' ')
    writetable(T3,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_stab3_nodamp.txt','Delimiter',' ')
    writetable(T4,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_stab4_nodamp.txt','Delimiter',' ')
    writetable(T5,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap6_SSIcov_on_stab5_nodamp.txt','Delimiter',' ')
end