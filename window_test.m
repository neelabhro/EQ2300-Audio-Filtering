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



    figure
    plot(nu(1:(N/2)), H_dB(1:(N/2)))
    hold on
    %plot(nu(1:(N/2)), H_quant_dB(1:(N/2)))
    ylim([-140 0])
    hold on
    
    
    
    
    
    
    
    
    
    
    
    
    
    
h_sinc = sinc(2*nu_c_low*(n-(M-1)/2));
A = 1/sum(h_sinc);
h_sinc = h_sinc * A;  % Normalize. H(0) = sum of h[n], so shoot for H(0) = 1

% Generate window
w = window(@bartlett, M);
%w = window(@hamming, M);
%w = window(@hanning, M);
%w = window(@blackman, M);
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

    plot(nu(1:(N/2)), H_dB(1:(N/2)))
    hold on
    %plot(nu(1:(N/2)), H_quant_dB(1:(N/2)))
    ylim([-140 0])
    hold on
    

    
    
    
    
    
    
    
    
h_sinc = sinc(2*nu_c_low*(n-(M-1)/2));
A = 1/sum(h_sinc);
h_sinc = h_sinc * A;  % Normalize. H(0) = sum of h[n], so shoot for H(0) = 1

% Generate window
%w = window(@bartlett, M);
%w = window(@hamming, M);
w = window(@hanning, M);
%w = window(@blackman, M);
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

    plot(nu(1:(N/2)), H_dB(1:(N/2)))
    hold on
    %plot(nu(1:(N/2)), H_quant_dB(1:(N/2)))
    ylim([-140 0])
    
    
        
h_sinc = sinc(2*nu_c_low*(n-(M-1)/2));
A = 1/sum(h_sinc);
h_sinc = h_sinc * A;  % Normalize. H(0) = sum of h[n], so shoot for H(0) = 1

% Generate window
%w = window(@bartlett, M);
%w = window(@hamming, M);
%w = window(@hanning, M);
%w = window(@blackman, M);
w = window(@chebwin, M);
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

    plot(nu(1:(N/2)), H_dB(1:(N/2)))
    hold on
    %plot(nu(1:(N/2)), H_quant_dB(1:(N/2)))
    xline(nu_c_low)
    xline(nu_c_high)
    yline(-40)
    ylim([-140 0])
    
    
    title('Magnitude of Frequency Response (dB)');
    xlabel('Normalized Frequency \nu');
    ylabel('20log_{10}|H_L(\nu)|');
    %legend('Blackman Original', 'Blackman quantization','Bartlett Original', 'Bartlett quantization','Hanning Original', 'Hanning quantization')    
    legend('Blackman','Bartlett','Hanning', 'Chebyshev') 
    %legend('Blackman','Hanning')  
    