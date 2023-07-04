% b_Stationarity    Calculate frame statistics and test stationarity of a
%                   random signal

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
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};

%--------------------------------------------------
% Start by running a frame statistics analysis using the RMS value
FileName='..\Data\RandomSignalExample';
load(FileName)
N=floor(length(y)/50);
F=framestat(y,N,'arms',0);
hf=figure;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 6]);
set(gcf, 'PaperSize', [12 6])
plot(F,'k','LineWidth',LineWidth)
xlabel('Frame number','FontName',FontName,'FontSize',FontSize)
ylabel('Running rms [m/s^2]','FontSize',FontSize)
A=axis;
axis([A(1) A(2) 0 1.5*A(4)])
grid
set(gca,'FontName',FontName,'FontSize',FontSize)

% Next, test the signal for stationarity with the reverse arrangements test
a=0.05;
status=teststat(F,a,'reverse');
if status
    fprintf('Data are stationary (reverse arrangements) with significance level %f\n',a)
else
    fprintf('Data are NOT stationary (reverse arrangements) with significance level %f\n',a)
end

% Next, test the signal for stationarity with the runs test
a=0.05;
status=teststat(F,a,'runs');
if status
    fprintf('Data are stationary (runs) with significance level %f\n',a)
else
    fprintf('Data are NOT stationary (runs) with significance level %f\n',a)
end
% To illustrate the difference between the reverse arrangements and runs
% test, we now add a periodic level modulation on the data, with one period
% during the length of the data, which is a nonstationary behavior.
% From Section 4.3 in Brandt, we know that the runs test should reveal
% this, but not the reverse arrangements test.
ttrend=linspace(0,2*pi,length(y));
ttrend=ttrend(:);                       % Force to column
F=framestat(y.*sin(ttrend),N,'std',0);
hf=figure;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 6]);
set(gcf, 'PaperSize', [12 6])
plot(F,'k','LineWidth',LineWidth)
xlabel('Frame number','FontName',FontName,'FontSize',FontSize)
ylabel('Running rms [m/s^2]','FontSize',FontSize)
A=axis;
axis([A(1) A(2) 0 1.5*A(4)])
grid
set(gca,'FontName',FontName,'FontSize',FontSize)

% Next, test the signal for stationarity with the reverse arrangements test
fprintf('********** Results from NONSTATIONARY data **********\n')
a=0.05;
status=teststat(F,a,'reverse');
if status
    fprintf('Data are stationary (reverse arrangements) with significance level %f\n',a)
else
    fprintf('Data are NOT stationary (reverse arrangements) with significance level %f\n',a)
end

% Next, test the signal for stationarity with the runs test
a=0.05;
status=teststat(F,a,'runs');
if status
    fprintf('Data are stationary (runs) with significance level %f\n',a)
else
    fprintf('Data are NOT stationary (runs) with significance level %f\n',a)
end

% Now, if you want to try one more thing on your own, add a trend to the
% data in y instead, and see that the reverse arrangements test then
% reveals the nonstationary behavior, but not the runs test.
% (Hint: the easiest trend is added to the variable F, by multiplying it
% with, for example, a linear trend from 1 to 2 during the length of F).
