function [f,H] = fft_fun(time,force,accln)

L = length(time);
% Sampling Freq
dt = time(2) - time(1);
Fs = 1/dt;
% Define a freq. vector
f = Fs/2*linspace(0,1,(L/2+1)); % in Hz
w = (f*2*pi)'; % in rad/sec
% Take the Fourier transform using fft
A = fft(accln,L);
A = A(1:length(accln)/2+1); % retain only the single sided spectrum
F = fft(force,L);
F = F(1:length(force)/2+1); % retain only the single sided spectrum
% Convert Accln to Disp
X = A./(-w.^2);
Sxx = conj(X).*X; % autopower output
Sff = conj(F).*F; % autopower output
Sxf = conj(X).*F; % crosspower
Sfx = conj(F).*X; % crosspower
% Calculate FRF
H = Sxf./Sff; % The correct way
H1 = X./F; % The quick way - may sometimes lead to errors
