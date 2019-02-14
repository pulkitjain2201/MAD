% Question 4
% Pulkit Jain
clc; clear; close all;

%% Part a(i)
load h1x1x_Endmill.txt
[F1,H11] = fft_fun(h1x1x_Endmill(:,1),h1x1x_Endmill(:,2),h1x1x_Endmill(:,3));
load h2x2x_Endmill.txt
[F2,H22] = fft_fun(h2x2x_Endmill(:,1),h2x2x_Endmill(:,2),h2x2x_Endmill(:,3));
load h3x3x_Endmill.txt
[F3,H33] = fft_fun(h3x3x_Endmill(:,1),h3x3x_Endmill(:,2),h3x3x_Endmill(:,3));
load h4x4x_Endmill_Medium.txt
[F4,H44] = fft_fun(h4x4x_Endmill_Medium(:,1),h4x4x_Endmill_Medium(:,2),h4x4x_Endmill_Medium(:,3));
