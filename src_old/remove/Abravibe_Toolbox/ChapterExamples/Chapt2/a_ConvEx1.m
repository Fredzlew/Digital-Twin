% a_ConvEx1         Illustration of convolution process


% This example produces the plot in Figure 2.9 of the book
% 

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

close all
clear

% Some settings to easily change plot appearance
%--------------------------------------------------
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
%--------------------------------------------------
% Set sampling frequency and create a time axis
fs=1000;
t=(0:1/fs:2000/fs)';

% Use the exponential decaying sine from Fig. 2.8, using the ABRAVIBE
% toolbox command aexpw (which is designed for impact testing, see Chapter
% 13)
e=aexpw(length(t),100*1e-4);        % End factor is in percent!
h=e.*sin(2*pi*2*t);
xmin=-3; xmax=3;
ymin=-1; ymax=1;
% Plot 
subplot(3,2,1)
plot(t,h,'-k','LineWidth',LineWidth)
hold on
plot([xmin xmax],[0 0],'-k','LineWidth',LineWidth)
plot([0 0],[ymin ymax],'-k','LineWidth',LineWidth)
xlabel('u','FontName',FontName,'FontSize',FontSize)
ylabel('h(u)','FontName',FontName,'FontSize',FontSize)
axis([xmin xmax ymin ymax])
set(gca,'FontName',FontName,'FontSize',FontSize)
set(gca,'XTick',[0],'FontName',FontName,'FontSize',FontSize)
set(gca,'XTickLabel',{'0'},'FontName',FontName,'FontSize',FontSize)
set(gca,'YTick',[0],'FontName',FontName,'FontSize',FontSize)
%grid
title('(a)','FontName',FontName,'FontSize',FontSize)
%-----------------------
% Flip the signal left-right, but it is in a column so this corresponds to
% flipping up/down
subplot(3,2,3)
hminusu=flipud(h);
t=t-2;
plot(t,hminusu,'-k','LineWidth',LineWidth)
hold on
plot([xmin xmax],[0 0],'-k','LineWidth',LineWidth)
plot([0 0],[ymin ymax],'-k','LineWidth',LineWidth)
xlabel('u','FontName',FontName,'FontSize',FontSize)
ylabel('h(-u)','FontName',FontName,'FontSize',FontSize)
axis([xmin xmax ymin ymax])
set(gca,'XTick',[0],'FontName',FontName,'FontSize',FontSize)
set(gca,'XTickLabel',{'0'},'FontName',FontName,'FontSize',FontSize)
set(gca,'YTick',[0],'FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
%grid
title('(b)','FontName',FontName,'FontSize',FontSize)
%-----------------------
% Shift to u=t, t=+1
subplot(3,2,5)
t=t+1;
htminusu=hminusu;
plot(t,htminusu,'-k','LineWidth',LineWidth)
hold on
plot([xmin xmax],[0 0],'-k','LineWidth',LineWidth)
plot([0 0],[ymin ymax],'-k','LineWidth',LineWidth)
xlabel('u','FontName',FontName,'FontSize',FontSize)
ylabel('h(t-u)','FontName',FontName,'FontSize',FontSize)
set(gca,'XTick',[0 1],'FontName',FontName,'FontSize',FontSize)
set(gca,'XTickLabel', {'0','t'},'FontName',FontName,'FontSize',FontSize)
set(gca,'YTick',[0],'FontName',FontName,'FontSize',FontSize)
axis([xmin xmax ymin ymax])
set(gca,'FontName',FontName,'FontSize',FontSize)
%text(1.1,-.3,'t','FontName',FontName,'FontSize',FontSize)
%grid
title('(c)','FontName',FontName,'FontSize',FontSize)
%-----------------------
% Right column: Convolution process
% Repeat lower left:subplot(3,2,5)
subplot(3,2,2)
plot(t,htminusu,'-k','LineWidth',LineWidth)
hold on
plot([xmin xmax],[0 0],'-k','LineWidth',LineWidth)
plot([0 0],[ymin ymax],'-k','LineWidth',LineWidth)
xlabel('u','FontName',FontName,'FontSize',FontSize)
ylabel('h(t-u)','FontName',FontName,'FontSize',FontSize)
set(gca,'XTick',[0 1],'FontName',FontName,'FontSize',FontSize)
set(gca,'XTickLabel', {'0','t'},'FontName',FontName,'FontSize',FontSize)
set(gca,'YTick',[0],'FontName',FontName,'FontSize',FontSize)
axis([xmin xmax ymin ymax])
set(gca,'FontName',FontName,'FontSize',FontSize)
%text(1.1,-.3,'t','FontName',FontName,'FontSize',FontSize)
%grid
title('(d)','FontName',FontName,'FontSize',FontSize)

%-----------------------
% e) x(u)
subplot(3,2,4)
tx=(0:1/fs:1)';
x=0.5*boxcar(length(tx));
x=[x;zeros(length(htminusu)-length(x),1)];
plot(t+1,x,'-k','LineWidth',LineWidth)
hold on
plot([xmin xmax],[0 0],'-k','LineWidth',LineWidth)
plot([0 0],[ymin ymax],'-k','LineWidth',LineWidth)
xlabel('u','FontName',FontName,'FontSize',FontSize)
ylabel('x(u)','FontName',FontName,'FontSize',FontSize)
set(gca,'XTick',[0 1],'FontName',FontName,'FontSize',FontSize)
set(gca,'XTickLabel', {'0','t'},'FontName',FontName,'FontSize',FontSize)
set(gca,'YTick',[0],'FontName',FontName,'FontSize',FontSize)
axis([xmin xmax ymin ymax])
set(gca,'FontName',FontName,'FontSize',FontSize)
%text(1.1,-.3,'t','FontName',FontName,'FontSize',FontSize)
%grid
title('(e)','FontName',FontName,'FontSize',FontSize)

%-----------------------
% f) area of x(u)h(t-u)
subplot(3,2,6)
Prod=2*flipud(x).*htminusu;
area(t,Prod)
%,'-k','LineWidth',LineWidth)
hold on
plot([xmin xmax],[0 0],'-k','LineWidth',LineWidth)
plot([0 0],[ymin ymax],'-k','LineWidth',LineWidth)
xlabel('u','FontName',FontName,'FontSize',FontSize)
ylabel('x(u) \cdot h(t-u)','FontName',FontName,'FontSize',FontSize)
set(gca,'XTick',[0 1],'FontName',FontName,'FontSize',FontSize)
set(gca,'XTickLabel', {'0','t'},'FontName',FontName,'FontSize',FontSize)
set(gca,'YTick',[0],'FontName',FontName,'FontSize',FontSize)
axis([xmin xmax ymin ymax])
set(gca,'FontName',FontName,'FontSize',FontSize)
%text(1.1,-.3,'t','FontName',FontName,'FontSize',FontSize)
%grid
title('(f) Area = y(t)','FontName',FontName,'FontSize',FontSize)
