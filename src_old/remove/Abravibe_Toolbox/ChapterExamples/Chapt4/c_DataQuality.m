% c_DataQuality       Data quality example
%
% Data quality example. Data in this example come from a test made on a
% truck. The 8 channels included, are some faulty, and some correct
% channels out of a larger set of data.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear
close all
clc

% Data directory
D='..\Data\TruckData\';
DirStruct=dir(strcat(D,'Truck*'));

% Define number of channels
NoChannels=length(DirStruct);

% First, run a frame statistics analysis of all channels
for n = 1:NoChannels
    FileName=strcat(D,DirStruct(n).name)
    load(FileName)
    N=floor(length(Data)/100);
    subplot(NoChannels/2,2,n)
    F(:,n)=framestat(Data,N,'std',1);
    A=axis;
    axis([A(1) A(2) 0 2*A(4)])
    title(['Frame Statistics, ''std'', Channel #' int2str(n)])
end

% As a comparison, we run a "reverse arrangements test" on channel 3
Frame=F(:,3);
NoFrames=100;
s=teststat(Frame,0.02,'reverse');
fprintf('******************************\n')
if s
    fprintf('Reverse arrangements test on channel 3 result: Data are stationary\n');
else
    fprintf('Reverse arrangements test on channel 3 result: Data are NOT stationary\n');
end
fprintf('******************************\n\n\n')

% Next, we run some standard statistics and log to a file
statchkf('..\Data\TruckData\Truck',1,8,'TruckStats');
% To list the file results
type TruckStats.log

% Finally, we plot the overall kurtosis of all channels, from the analysis just
% made
load TruckStats
Kurtosis=Kurtosis/Kurtosis(3);           % Normalize kurtosis to 'channel' 3
figure
bar(Kurtosis,'FaceColor',[0.8 0.8 0.8])
xlabel('File Number')
title('Kurtosis normalized by ch. #3')
grid

