%% Script that investigate each sensor %%
clear;clc;close all
data = readmatrix('data_5_6_1.txt')'; % Loading displacement data
fss = data(2:6,:)/1000; % Converting mm to m
xm = flip(fss,1); % Swap rows due to sensor

xm_std = zeros(5,1);
xm_mean = zeros(5,1);

for i = 1:5
xm_std(i) = std(xm(i,:));
xm_mean(i) = mean((xm(i,:)));
end

xbase = data(7,:)/1000; % Converting mm to m
xm = flip(fss,1)-xbase; % Swap rows due to sensor
