%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EQ 2300 - Digital Signal Processing
% All tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Select whether or not to show plots
plot_windows = 1;
plot_impulse_responses = 1;
plot_frequency_responses = 1;

nu_c_low  = 1/16;
nu_c_high = 1/8;

N = 1024;       % Samples for FFT. If higher than M, will pad with zeros. The higher N, the finer sampling of H(nu)
M = 47;         % Number of taps in the filter. MUST BE ODD in order to get a Type I FIR filter
F = 14;         % Number of bits used to quantize FIR coefficients

n = 0:1:M-1;
nu = linspace(0,1,N);

% Generate sinc impulse response (ideal LPF)
h_sinc = sinc(2*nu_c_low*(n-(M-1)/2));
A = 1/sum(h_sinc);
h_sinc = h_sinc * A;  % Normalize. H(0) = sum of h[n], so shoot for H(0) = 1

% Generate window
%w = window(@bartlett, M);
%w = window(@hamming, M);
%w = window(@hanning, M);
w = window(@blackman, M);
%w = window(@chebwin, M);
%w = window(@kaiser, M);

% Trucante ideal LPF with window
h = zeros(1,M);
for i = 1:M
    h(i) = h_sinc(i) * w(i);
end
H = fft(h,N);
H_dB = 20 * log10(abs(H));

% High-pass filter
delta = zeros(1,M);
delta((M+1)/2) = 1; % delta has to be delayed (M+1)/2 samples so h_HPF is symmetric, and therefore a valid type I FIR filter
h_HPF = delta - h;
H_HPF = fft(h_HPF,N);
H_HPF_dB = 20 * log10(abs(H_HPF));


% Apply quantization effect
h_quant = 2^(-F) * round(h * 2^F);
H_quant = fft(h_quant, N);
H_quant_dB = 20 * log10(abs(H_quant));




% Task 5 - SQNR
% Power of the signal after h_Q[n], still without the effect of quantization
P_x = 2^22/3;
P_xL = P_x * sum(h.^2);
% Power of quantization noise
step = 2^11 / 2^(F-1);
P_q = step^2 / 12;
SQNR = P_xL / P_q;
SQNR_dB = 10 * log10(SQNR)

if plot_windows    
    figure
    stem(n, w, 'filled')
    title('Window used')
end

if plot_impulse_responses
    figure
    stem(n, h, 'filled')
    title('Sinc truncated with window')
    figure
    stem(n,h_HPF, 'filled')
    title('HPF')
end

if plot_frequency_responses
    figure
    plot(nu(1:(N/2)), H_dB(1:(N/2)))
    hold on
    plot(nu(1:(N/2)), H_quant_dB(1:(N/2)))
    xline(nu_c_low, 'red')
    xline(nu_c_high, 'red')
    yline(-40, 'red')
    ylim([-140 0])
    legend('Original response', 'After quantization')
    
    figure
    plot(nu(1:(N/2)), H_HPF_dB(1:(N/2)))
    title('HPF')
    xline(nu_c_low, 'red')
    xline(nu_c_high, 'red')
    ylim([-60 0])
end
