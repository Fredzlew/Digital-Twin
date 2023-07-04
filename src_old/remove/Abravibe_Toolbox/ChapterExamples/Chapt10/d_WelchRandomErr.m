% d_WelchRandomErr example of random error in a PSD estimated using Welch's 
%                  method using different overlap percentages.
%

% This example is similar to the plot in Fig. 10.8 in Brandt.

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
% Relative increase of equivalent number averages
Xaxis=[0:10:90];
N=8192;
Nsamp=100*N;
w=[boxcar(N) ahann(N) bartlett(N)];
for wn = 1:length(w(1,:))
    for k = 1:length(Xaxis)
        blkstep=N-floor(Xaxis(k)*N/100);
        K(k)=floor((Nsamp-N)/blkstep+1);
        rerrt(k,wn)=welcherr(w(:,wn),Xaxis(k),K(k));    % Theoretical rand err
    end
end
Me=1./rerrt.^2;
plot(Xaxis,Me(:,2)/100,LineType{1},Xaxis,Me(:,3)/100,LineType{2},Xaxis,Me(:,1)/100,LineType{3},'LineWidth',LineWidth)
xlabel('Overlap [%]','FontName',FontName,'FontSize',FontSize)
ylabel('Relative increase of equivalent number averages','FontName',FontName,'FontSize',FontSize)
grid
axis([0 90 1 2.5])
legend('Hanning','Bartlett','Rectangular')
set(gca,'FontName',FontName,'FontSize',FontSize)
