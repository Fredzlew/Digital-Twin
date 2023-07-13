% parameters
clc; clear; close all;
addpath(genpath('..\data'),genpath('..\npy-matlab-master'))
SSIstab2 = readNPY('..\data\simulated_data\Modal_par\SSI_stab_sim.npy');

% Only use frequency up to 12 
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

T1 = array2table([num2cell(freq0'),num2cell((modelorder0*2)'),num2cell(damp0')]);
T2 = array2table([num2cell(freq1'),num2cell((modelorder1*2)'),num2cell(damp1')]);
T3 = array2table([num2cell(freq2'),num2cell((modelorder2*2)'),num2cell(damp2')]);
T4 = array2table([num2cell(freq3'),num2cell((modelorder3*2)'),num2cell(damp3')]);
T5 = array2table([num2cell(freq4'),num2cell((modelorder4*2)'),num2cell(damp4')]);

T1.Properties.VariableNames(1:3) = {'freq0','modelorder0','damp0'};
T2.Properties.VariableNames(1:3) = {'freq1','modelorder1','damp1'};
T3.Properties.VariableNames(1:3) = {'freq2','modelorder2','damp2'};
T4.Properties.VariableNames(1:3) = {'freq3','modelorder3','damp3'};
T5.Properties.VariableNames(1:3) = {'freq4','modelorder4','damp4'};

writetable(T1,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap4_SSIcov_on_stab1_sim.txt','Delimiter',' ')
writetable(T2,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap4_SSIcov_on_stab2_sim.txt','Delimiter',' ')
writetable(T3,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap4_SSIcov_on_stab3_sim.txt','Delimiter',' ')
writetable(T4,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap4_SSIcov_on_stab4_sim.txt','Delimiter',' ')
writetable(T5,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap4_SSIcov_on_stab5_sim.txt','Delimiter',' ')