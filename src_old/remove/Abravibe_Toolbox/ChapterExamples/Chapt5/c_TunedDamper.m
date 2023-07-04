% c_TunedDamper     Compute FRF on SDOF system with tuned damper added
%
% This example simulates addition of a tuned damper on a SDOF system. It
% uses the ABRAVIBE command tunedamp to simulate the FRF before and after
% adding the damper.


clear
clc
close all

% Some settings to easily change plot appearance
%--------------------------------------------------
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
LineType={'-','--','-.',':'};
%--------------------------------------------------

[H,f]=tunedamp(100,.003,.1,10);
semilogy(f,abs(H(:,1)),'-',f,abs(H(:,2)),'--','LineWidth',LineWidth);
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Dynamic Flexibility [m/N]','FontName',FontName,'FontSize',FontSize)
grid
set(gca,'FontName',FontName,'FontSize',FontSize)
legend('FRF w. tuned damper','FRF of SDOF only')
