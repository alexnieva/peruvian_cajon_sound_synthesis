% Coupled Mode Synthesis for the Cajon drum by Alex Nieva.
% http://alexnieva.github.io/
% Final project MUMT 618 Fall 2015 McGill University.
% Prof. Gary Scavone.

clear all; close all;


% Sample or real cajon damped - low frequency hit
% Change to 'CAJON_SAMPLE3.wav' for high frequency hit
% These measurements were done strumming the cajon with the hand in 
% the way a regular cajon player will do it.
% Unfortunately, the cajon used had snares attached to the back of the
% front plates and only sound pressure measurements could be performed.
[x, fs]=audioread('CAJON_SAMPLE1.wav');         
%[xx,fs]=audioread('CAJON_SAMPLE_BACK_UD1.wav'); % Sample or real cajon with snares 
                                                % recorded in output hole
[SB,fs]=audioread('CAJON_SAMPLE_SUSANABACA3low.wav'); % Sample peruvian cajon 
N = length(x);
ximp = [1 zeros(1, N-1)]; % Unit impuse for testing
%x = ximp; % Uncomment this line to get impulse response

%fs = 44100;
T = 1 / fs;

% first reson
f0 = 70; % Resonant frequency
B = 4; % Bandwidth Hertz
Brad = 2*pi*B/fs;
R = 1 - Brad/2; %pole of the filter.
f0rad = 2*pi*f0/fs;
cosine_theta = cos(f0rad)*(2*R)/(1+R^2);
theta = acos(cosine_theta);
A0 = (1-R^2)*sin(theta);

% second reson
f1 = 100; % Resonant frequency
B1 = 10; % Bandwidth Hertz
B1rad = 2*pi*B1/fs;
R1 = 1 - B1rad/2;
f1rad = 2*pi*f1/fs;
cosine_theta1 = cos(f1rad)*(2*R1)/(1+R1^2);
theta1 = acos(cosine_theta1);
A1 = (1-R1^2)*sin(theta1);

% Third reson (second harmonic of helmholtz reson)
f2 = 180; % Resonant frequency
B2 = 10; % Bandwidth Hertz
B2rad = 2*pi*B2/fs;
R2 = 1 - B2rad/2;
f2rad = 2*pi*f2/fs;
cosine_theta2 = cos(f2rad)*(2*R2)/(1+R2^2);
theta2 = acos(cosine_theta2);
A2 = (1-R2^2)*sin(theta2);

% fourth reson (plate mode - approximate)
f3 = 500; %Resonant frequency
B3 = 80; %Hertz
B3rad = 2*pi*B3/fs;
R3 = 1 - B3rad/2; %pole of the filter.
f3rad = 2*pi*f3/fs;
cosine_theta3 = cos(f3rad)*(2*R3)/(1+R3^2);
theta3 = acos(cosine_theta3);
A3 = (1-R3^2)*sin(theta3);

% fifth reson (for higher frequency content)
f4 = 2050; %Resonant frequency
B4 = 1000; %Hertz
B4rad = 2*pi*B4/fs;
R4 = 1 - B4rad/2; %pole of the filter.
f4rad = 2*pi*f4/fs;
cosine_theta4 = cos(f4rad)*(2*R4)/(1+R4^2);
theta4 = acos(cosine_theta4);
A4 = (1-R4^2)*sin(theta4);

% Filters
A = [A0 0 -1]; % Added zeros to force low frequency decay.
B = [1 -2*R*cosine_theta R^2];
C = [A1 0 -1]; % Added zeros to force low frequency decay.
D = [1 -2*R1*cosine_theta1 R1^2];
E = [A2 0 -1];
F = [1 -2*R2*cosine_theta2 R2^2];
G = A3;
H = [1 -2*R3*cosine_theta3 R3^2];
I = A4;
J = [1 -2*R4*cosine_theta4 R4^2];

sound(x,fs);
disp('Sound of damped hit');
pause;

y = filter(A,B,x);
y = 0.95 * y / max(abs(y));
figure; subplot(311); plot(y)
sound(y, fs)
disp('Synthesized sound of box as Helmholtz Resonator');
pause;

z = filter(C,D,x);
z = 0.35 * z / max(abs(z));
subplot(312); plot(z,'r');

p = filter(E,F,x);
p = 0.25 * p / max(abs(p));
subplot(313); plot(p);
sound((y+z+p)./(max(abs(y+z+p))), fs)
disp('Synthesized sound of the box')
pause;

plate1 = filter(G,H,x);
plate1 = 0.5*plate1/max(abs(plate1));

plate2 = filter(I,J,x);
plate2 = 0.6*plate2/max(abs(plate2));

% Playing sounds
output = (y+z+p+plate1+plate2)./(max(abs(y+z+p+plate1+plate2)));
sound(output, fs)
disp('Synthesized sound of Cajon including plate resonances')
pause;

sound(SB,fs);
disp('Sound of real sampled Cajon - two hits (taken from audio recorded file)')

% Fourier Transform Analysis to see frequency content.
Nfft = 2^14;
ft1 = fft(x,Nfft);
ft2 = fft(output,Nfft);
ft3 = fft(SB,Nfft);

ff = 0:fs/Nfft:(fs - fs/Nfft);

% % Impulse response UNCOMMENT WHEN USING ximp!!!
% IR = fft(output,Nfft);
% figure; semilogx(ff,20*log10(abs(IR)));
% axis([10 2500 0 70]); grid
% title('Frequency Response Synthesis Filters');
% xlabel('frequency(Hz)')
% ylabel('dB')
% 

figure;
subplot(311);plot(ff,abs(ft1)); axis([0 2500 0 2000])
title('Ft of input (dry sound)')
subplot(312);plot(ff,abs(ft2)); axis([0 2500 0 2000])
title('FT Synthesized sound')
subplot(313);plot(ff,abs(ft3)); axis([0 2500 0 2000])
title('FT Real Cajon (two hits)')

