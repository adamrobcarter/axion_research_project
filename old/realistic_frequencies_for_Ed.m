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

% plot directly
%{
figure
subplot(2, 1, 1)
plot(t, f_uncorrelated)
subplot(2, 1, 2)
plot(t, f_correlated)
uiwait
%}

% plot histograms
%{
figure
subplot(2, 1, 1)
histogram(f_uncorrelated)
subplot(2, 1, 2)
histogram(f_correlated)
uiwait
%}

% the moving mean averaging clearly reduces the std. deviation of the
% points, but by how much? must be a theoretical answer, or could just
% calculate from the data. By inspection, it appears to be reduced by a 
% factor of sqrt(boxcar_length)
std_dev_uncorrelated = std(f_uncorrelated);
std_dev_correlated = std(f_correlated);

% plot autocorrelation
%{
figure
subplot(2, 1, 1)
[acf, k] = autocorr(f, NumLags=num_timesteps-1);
plot(k/f_sampling, acf)
subplot(2, 1, 2)
[acf, k] = autocorr(f_correlated, NumLags=num_timesteps-1);
plot(k/f_sampling, acf)
uiwait
%}

%% turn into an array of frequencies
% renormalise f_correlated to have std.dev = 1
f_renormalised = f_correlated / std_dev_correlated;
% confirm this worked
std_dev_renormalised = std(f_renormalised);
% now do the axion stuff, commented out while we work just on noise

f_0 = 0.1* axion_rest_energy / h_eV; % Hz, default axion frequency
axion_linewidth = f_0 * 0.5 * (axion_velocity / c)^2;

% array in time of correlated frequencies with given linewidth
f_axion = f_0 + axion_linewidth * f_renormalised;

f_axion = linspace(f_0/10, f_sampling/2, length(t));

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

%SNR = power_signal / power_noise * sqrt( axion_linewidth * t_end);

%{
%% run iwave
% calculate avg df/dt (for calculating tau_iwave)
dfdsample = f(1:end-1) - f(2:end);
dfdt = dfdsample * f_sampling;
dfdt_square = dfdt .* dfdt;
dfdt_square_avg = mean(dfdt_square);

tau_iwave = (288 * pi^4 * f_0^2 * dfdt_square_avg )^(-1/6); % tau_opt from the iwave paper
tau_iwave = 1.5e-6;

f_iwave = extract_frequency(f_sampling, f_0, t, x, tau_iwave);

% plot iwave before and after

figure
plot(t, f)
lim = ylim();
hold on
plot(t, f_iwave)
ylim(lim)
legend('before iwave', 'after iwave')
xlabel('time (s)')
ylabel('frequenzy (Hz)')
uiwait
%}

%% now feed our wave into a resonator
x = signal;
x = noise; % skip wave for the mo and just act on noise
Q = 5e3;

f = f_0 * ones(size(t));
y = resonator(x, f_0, f, Q, f_sampling);
f = display_spectra(x, y, f_sampling);
hold on

omega_0 = 2*pi * f_0;
Delta_0 = omega_0 / f_sampling;
w = Delta_0 / (2 * Q);
Delta = 2*pi*f/f_sampling;
A_Delta = (1 - exp(-w)) ./ sqrt(1 - 2*exp(-w)*cos(Delta - Delta_0) + exp(-2*w)); % from Ed paper
plot(f, A_Delta, 'DisplayName','A(Delta)')

s = 1i * 2 * pi * f;
Gamma = omega_0 / Q;
Hs = Gamma * s ./ ( s.^2 + Gamma*s + omega_0^2 ); % also from Ed paper
plot(f, abs(Hs), 'DisplayName','|H(s)|')

uiwait

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
function f = display_spectra(x, y, f_sampling)
    [amp_bef, f] = nsd_pwelch(x, 0.000001, f_sampling);
    [amp_aft, f] = nsd_pwelch(y, 0.000001, f_sampling);

    subplot(1,2,1)
    loglog(f, amp_bef)
    hold on 
    loglog(f, amp_aft)
    legend('before', 'after', 'Location', 'northwest')
    xlim([1e7, f_sampling/2])
    title('simulated power spectrum before/after resonant filter')

    subplot(1,2,2)
    ratio = amp_aft./amp_bef;
    loglog(f, ratio)
    legend('simulation', 'Location', 'northwest')
    title('filter transfer function')
    xlim([1e7, f_sampling/2])
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
function freq = extract_frequency(f_sampling, f_guess, time, wave, tau)
    % initialise iwave pll
    start_time = 0.0; % secs
    fstate = iwavePllConstruct(f_sampling, tau, f_guess, start_time);
    
    % run iwave pll
    [dout, qout, amp, error, freq, fstate] = iwavePllRun(wave, fstate);
end