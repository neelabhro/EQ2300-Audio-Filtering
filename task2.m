%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EQ 2300 - Digital Signal Processing
% Task 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close;
clc;

nu_c_low  = 1/16;
nu_c_high = 1/8;

N = 100;
M = 25;
n = 0:1:N-1;

% Generate sinc impulse response
h = sinc(2*nu_c_low*(n-M));
A = 1/sum(h);
h = h * A;  % Normalize. H(0) = sum of h[n], so shoot for H(0) = 1
stem(h);
xlabel('n');
ylabel('h[n]');
figure;

% Ideal LPF
N_fft = 1024;
nu = linspace(0,1,N_fft);
H_ideal = [ones(1, round(N_fft*nu_c_low)), zeros(1, N_fft - round(N_fft*nu_c_low))];
H_ideal_dB = 20 * log10(H_ideal);

% Realizable LPF, from truncated sinc
H = fft(h,N_fft); % FFT of a normalized sinc impulse response
H_dB = 10 * log10(abs(H));

plot(nu(1:N_fft/2), abs(H(1:N_fft/2)));
hold on
plot(nu(1:N_fft/2), H_ideal(1:N_fft/2));
title('Impulse Response');
xlabel('n');
ylabel('h_L[n]');
figure;

plot(nu(1:(N_fft/2)), H_dB(1:(N_fft/2)));
hold on
plot(nu(1:(N_fft/2)), H_ideal_dB(1:(N_fft/2)));
title('Magnitude of Frequency Response');
xlabel('Normalized Frequency nu');
ylabel('20log_{10}|H_L(nu)|');
figure;



%% Plot filter freq response - compare with design specifications

plot(nu(1:(N_fft/2)), H_dB(1:(N_fft/2)))
xline(nu_c_low, 'red')
xline(nu_c_high, 'red')
yline(-40, 'red')
title('Filter Frequency Response');
figure;

%% Windows

M = 51;

w_bartlett = window(@bartlett,M);
w_hamming = window(@hamming, M);
w_chebysev = window(@chebwin, M);

subplot(1,3,1)
stem(w_bartlett, 'filled')
title('Bartlett')
subplot(1,3,2)
stem(w_hamming, 'filled')
title('Hamming')
subplot(1,3,3)
stem(w_chebysev, 'filled')
title('Chebysev')
set(gcf,'Position',[0 300 1200 300])

figure;
%% Compare windows
% Bartlett's window

h_w = zeros(1,M);
for i = 1:M
    h_w(i) = h(i) * w_bartlett(i);
end
%stem(h_w,'filled')

N_fft = 1024;
nu = linspace(0,1,N_fft);
H_w = fft(h_w,N_fft);
H_w_dB = 20 * log10(abs(H_w));


subplot(1,3,1)
plot(nu(1:(N_fft/2)), H_w_dB(1:(N_fft/2)))
title('Bartlett')
xline(nu_c_low, 'red')
xline(nu_c_high, 'red')
yline(-40, 'red')
ylim([-140 0])

% Hamming's window

h_w = zeros(1,M);
for i = 1:M
    h_w(i) = h(i) * w_hamming(i);
end
%stem(h_w,'filled')

N_fft = 1024;
nu = linspace(0,1,N_fft);
H_w = fft(h_w,N_fft);
H_w_dB = 20 * log10(abs(H_w));

subplot(1,3,2)
plot(nu(1:(N_fft/2)), H_w_dB(1:(N_fft/2)))
title('Hamming')
xline(nu_c_low, 'red')
xline(nu_c_high, 'red')
yline(-40, 'red')
ylim([-140 0])

% Chebysevs's window

h_w = zeros(1,M);
for i = 1:M
    h_w(i) = h(i) * w_chebysev(i);
end
%stem(h_w,'filled')

N_fft = 1024;
nu = linspace(0,1,N_fft);
H_w = fft(h_w,N_fft);
H_w_dB = 20 * log10(abs(H_w));

subplot(1,3,3)
plot(nu(1:(N_fft/2)), H_w_dB(1:(N_fft/2)))
title('Chebysev')
xline(nu_c_low, 'red')
xline(nu_c_high, 'red')
yline(-40, 'red')
ylim([-140 0])
set(gcf,'Position',[0 300 1200 300])