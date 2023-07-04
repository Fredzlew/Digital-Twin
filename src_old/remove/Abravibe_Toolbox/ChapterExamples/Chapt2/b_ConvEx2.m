% b_ConvEx2         Illustration of the convolution process

% This example produces a plot similar to Figure 2.10 in the book

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
% Load some example data from a truck, and create a time axis
% File includes variables y and fs
load ..\Data\RandomSignalExample
tx=makexaxis(y,1/fs);
x=y;                    % Rename signal to input x(t)
% Create a time axis for a 0.3 seconds impulse response
th=(0:1/fs:.3)';
% Use an exponential decaying sine 
e=aexpw(length(th),1);        % End factor is in percent
h=0.02*e.*sin(2*pi*18*th);    % Arbitrary scaling to get result within reasonable amplitude
% Use h(t) to convolve up to three seconds. x is approx. 5 seconds long
tidx=find(tx<=3);
y=conv(x,h);
y=y(tidx);                   % Conv produces a little longer signal; truncate it, 
                             % as we want to illustrate the principle only
ty=tx(tidx);
xmin=min(tx); xmax=max(tx);
ymin=-5; ymax=5;
% Plot 
subplot(3,1,1)
N=length(h);
th=(max(ty):-1/fs:max(ty)-(N-1)/fs)';       % 'Backwards' x axis for h!
plot(th,h,'-k','LineWidth',LineWidth)
axis([xmin xmax min(h) max(h)])
ylabel('h(3-u)','FontName',FontName,'FontSize',FontSize)
set(gca,'YTick',[])
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(3,1,2)
hold on
box('on');
area([min(th) max(th)],[5 5],'FaceColor',[0.8 0.8 0.8],'BaseValue',-5);
plot(tx,x,'-k','LineWidth',LineWidth)
axis([xmin xmax ymin ymax])
ylabel('x(t)','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(3,1,3)
plot(ty,y,'-k','LineWidth',LineWidth)
hold on
plot(ty(end),y(end),'*k','LineWidth',LineWidth)
axis([xmin xmax ymin ymax])
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('y(t)','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
