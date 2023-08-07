clc; clear; close all;

FDDPSD = readNPY('..\data\simulated_data\Modal_par\FDDPSD.npy');
FDDPSDfreq = readNPY('..\data\simulated_data\Modal_par\FDDPSDfreq.npy');

% Taking the PSD out of the 3D matrix
for i = 1:length(FDDPSD)
    PSD1 = diag(FDDPSD(:,:,i));
    PSD(i,1) = 10*log(PSD1(1));
    PSD(i,2) = 10*log(PSD1(2));
    PSD(i,3) = 10*log(PSD1(3));
    PSD(i,4) = 10*log(PSD1(4));
    PSD(i,5) = 10*log(PSD1(5));
end

% load frequencies
f_ = FDDPSDfreq(1:length(FDDPSD));

figure
for i = 1:5
hold on
plot(f_,(PSD(:,i)))
end
xlim([0 12])
hold off

%% Download the file to plot in latex
%{
T = array2table([num2cell(f_(1:1201)),num2cell(PSD(1:1201,:))]);
T.Properties.VariableNames(1:6) = {'f_','PSD1','PSD2','PSD3','PSD4','PSD5'};
writetable(T,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap4_FDD_PSD_sim.csv','Delimiter',';')
%}