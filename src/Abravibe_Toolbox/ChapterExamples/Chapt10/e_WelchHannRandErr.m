% e_WelchHannRandErr    Plot random error as function of number of averages,
%                       when using Welch with Hanning and 50% overlap.
%

% This example is similar to the plot in Fig. 10.9 in Brandt.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear
clc
close all

% Some settings to easily change plot appearance
%--------------------------------------------------
FontSize=11;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};

%--------------------------------------------------
% Calculate normalized random error for different number of total FFTs
NBlocks=[1e1 1e2 1e3 1e4];
for n=1:length(NBlocks)
    er(n)=welcherr(ahann(8192),50,NBlocks(n));
end
% Plot results
loglog(NBlocks,er,LineType{1},'LineWidth',LineWidth)
xlabel('Number of FFTs','FontName',FontName,'FontSize',FontSize)
ylabel('Normalized random error,    \epsilon_r','FontName',FontName,'FontSize',FontSize)
grid
axis([10 1e4 1e-2 .5])
set(gca,'FontName',FontName,'FontSize',FontSize)
