f_sampling = 2e9; % Hz
t_end = 0.01; % secs

t_step = 1/f_sampling;
t = 0 : t_step : t_end;

%% constants
c = 3e8; % m/s
hbar_times_c = 2e-7; % eV m
h_eV = 4.136e-15; % eV
k_B = 1.38065e-23; % m2 kg s-2 K-1

%% parameters
axion_rest_energy = 2e-6; % eV/c2
axion_velocity = 220e3; % m/s
T = 10e-3; % K
axion_power = 1.52e-21; % W, from eq 2.28 of Daw thesis

%% uncorrelated gaussian points
f_uncorrelated = randn(size(t));

%% correlate the points
de_broglie_wavelength = 2*pi * hbar_times_c / (axion_rest_energy * axion_velocity/c);
t_coherence = de_broglie_wavelength / axion_velocity;

correlation_time = t_coherence;
assert(correlation_time < t_end / 2) % this not being true caused a long problem once
boxcar_length = correlation_time * f_sampling; % samples
f_correlated = movmean(f_uncorrelated, boxcar_length);
clear f_uncorrelated

%% turn into an array of frequencies
% renormalise f_correlated to have std.dev = 1
f_correlated = f_correlated / std(f_correlated);
% confirm this worked
std_dev_renormalised = std(f_correlated);
% now do the axion stuff, commented out while we work just on noise

f_0 = axion_rest_energy / h_eV; % Hz, default axion frequency
axion_linewidth = f_0 * 0.5 * (axion_velocity / c)^2;
Q_axion = f_0 / axion_linewidth;

% array in time of correlated frequencies with given linewidth
f_axion = f_0 + axion_linewidth * f_correlated;
clear f_correlated

%% calculate a wave in time with the given frequency
signal = calc_wave_with_frequency(t, t_step, f_axion);
signal_power = bandpower(signal); % = 0.5000 - must be a reason for this exactness
signal = signal / sqrt(signal_power) * sqrt(axion_power);
power_signal = bandpower(signal);

%% calculate thermal noise;
noise = randn(size(t));
% this noise has total power 1, now scale it
noise = noise * sqrt(k_B * T * f_sampling / 2);
power_noise = bandpower(noise);

%% now feed our wave into a resonator

f = f_0 * ones(size(t));

Q = [1e4 1e5 1e6 1e7 1e8 1e9 1e10];
%a = [0 0 0 0 0 0 0];

for i = 1:7
    %a(i) = get_static_SNR(signal, noise, f_0, f, Q(i), f_sampling)
end

scatter(Q/Q_axion, a)
set(gca,'xscale','log','yscale','log')
xlabel('Q/Q_a')
ylabel('power ratio')
adstyle(8, 8)
saveas(gcf,'figures/static_q_limit.jpg')

function static_SNR = get_static_SNR(signal, noise, f_0, f, Q, f_sampling)
    static_signal = resonator(signal, f_0, f, Q, f_sampling);
    static_noise  = resonator(noise,  f_0, f, Q, f_sampling);
    static_SNR = bandpower(static_signal) / bandpower(static_noise);
end

%% end

% pass the supplied wave x through a resonator
% given centre frequency f_0_centre
% resonant frequency (array same size as x) f_0,
% Q, and f_sampling
function y = resonator(x, f_0_centre, f_0, Q, f_sampling)
    % this all from Daw 2019 appendix
    % set up params
    omega_0 = 2*pi * f_0;
    Delta_0 = omega_0 / f_sampling;
    w = Delta_0 / (2 * Q);  
    
    % find the quadratures
    I = x;
    quarter_period_secs = 1/4 * 1/f_0_centre;
    quarter_period_samples = quarter_period_secs * f_sampling;
    Q = delay(x, quarter_period_samples);

    x_q = I + 1i * Q;
    y_q = zeros(size(x_q));
    
    % do the algorithm
    for i = 2 : length(x_q) % should it be Delta_0(i) or (i-1)?
        y_q(i) = exp(-w(i) + 1i * Delta_0(i))*y_q(i-1) + (1 - exp(-w(i)))*x_q(i);
    end
    
    % recombine the quadratures
    I = delay(real(y_q), quarter_period_samples);
    Q = imag(y_q);
    
    y = I + Q;
end

% delay a signal by the specified number of samples.
% pads the start with zeroes
function x_delayed = delay(x, delay_samples)
    delay_samples = floor(delay_samples); % needs to be integer
    x_delayed = zeros(size(x));
    x_delayed(delay_samples+1:end) = x(1:end-delay_samples);
end

% show two spectra and their transfer function
function [f, ratio] = display_spectra(x, y, f_sampling)
    [amp_bef, f] = nsd_pwelch(x, 0.000005, f_sampling);
    [amp_aft, f] = nsd_pwelch(y, 0.000005, f_sampling);

    subplot(1,2,1)
    loglog(f, amp_bef)
    hold on 
    loglog(f, amp_aft)
    legend('before', 'after', 'Location', 'southwest')
    xlim([1e7, f_sampling/2])
    title('simulated power spectrum before/after resonant filter')

    subplot(1,2,2)
    ratio = amp_aft./amp_bef;
    loglog(f, ratio)
    legend('simulation', 'Location', 'southwest')
    title('filter transfer function')
    xlim([1e7, f_sampling/2])

    [pks,locs,w,p] = findpeaks(ratio, WidthReference="halfheight", SortStr='descend');
    maxindex = 1;
    peak_freq = f(locs(maxindex));
    peak_w = w(maxindex)
    %peak_Q = peak_freq / peak_w
end

% calculate a wave in time with the given frequency in time
function x = calc_wave_with_frequency(t, t_step, f)
    x = zeros(size(t));
    f_integral = 0;

    for i = 1 : length(t)
        x(i) = sin(2*pi*f_integral);
        f_integral = f_integral + f(i)*t_step;
    end
end

% extract the frequency from a wave using IWAVE
function freq = extract_frequency(wave, f_sampling, f_guess, tau)
    % initialise iwave pll
    start_time = 0.0; % secs
    fstate = iwavePllConstruct(f_sampling, tau, f_guess, start_time);
    
    % run iwave pll
    [dout, qout, amp, error, freq, fstate] = iwavePllRun(wave, fstate);
end