function improvement = realistic_frequencies()
f_sampling = 2e9; % Hz
t_end = 0.01; % secs

t_step = 1/f_sampling;
t = 0 : t_step : t_end;

%% constants
c = 3e8; % m/s
h_eV = 4.136e-15; % eV
k_B = 1.38065e-23; % m2 kg s-2 K-1
h = 6.62607015e-34; % m2 kg / s
hbar_times_c = 2e-7; % eV m

%% parameters
axion_rest_energy = 2e-6; % eV/c2
axion_velocity = 220e3; % m/s
T = 10e-3; % K
axion_power = 1.52e-21; % W, from eq 2.28 of Daw thesis

%% axion signal params
f_0 = axion_rest_energy / h_eV; % Hz, default axion frequency
axion_linewidth = f_0 * 0.5 * (axion_velocity / c)^2;
axion_Q = f_0 / axion_linewidth;
de_broglie_wavelength = 2*pi * hbar_times_c / (axion_rest_energy * axion_velocity/c);
t_coherence = de_broglie_wavelength / axion_velocity;

%% generate signal and noise
signal = generate_axion_signal(t, f_sampling, f_0, axion_power, axion_linewidth, t_coherence);
noise = generate_thermal_noise(t, f_sampling, T);

%% run iwave
tau_iwave = 1e-5; % before adding noise
tau_iwave = 3e-4;
tau_iwave = 5e-4; % for 1* mult
%tau_iwave = 20e-4; % for 10* mult

% atm we're running iwave on the signal not the signal + noise
iwave_input = signal + noise;

band_width = 150;
band = [f_0-band_width/2, f_0+band_width/2];
band = [f_0-band_width*0.2, f_0+band_width*0.8];
iwave_input = bandpass(iwave_input, band, f_sampling);

noise_power_in_band = bandpower(noise, f_sampling, band);
signal_power_in_band = bandpower(signal, f_sampling, band);
combined_power_in_band = bandpower(iwave_input, f_sampling, band);
combined_power_total = bandpower(iwave_input);

f_iwave = run_iwave(iwave_input, f_sampling, f_0, tau_iwave);

% now remove the start as it can go a bit off
start_secs = 5 * tau_iwave;
start_secs = 0.000001;
start_samples = start_secs * f_sampling;
t       = t      (start_samples:end);
f_iwave = f_iwave(start_samples:end);
f_axion = f_axion(start_samples:end);
signal  = signal (start_samples:end);
noise   = noise  (start_samples:end);

% find avg diff between iwave and original
avg_dist_dynamic = rms(f_iwave - f_axion);
Qmax_dynamic = 2 * f_0 / avg_dist_dynamic;

avg_dist_static = rms(f_0 * ones(size(t)) - f_axion);
Qmax_static = 2 * f_0 / avg_dist_static;

% plot iwave before and after

figure
plot(t, f_axion)
lim = ylim();
hold on
plot(t, f_iwave)
ylim(lim)
legend('before iwave', 'after iwave')
xlabel('time (s)')
ylabel('frequenzy (Hz)')
xline(t_coherence)
xline(tau_iwave)

%% now feed our wave into a resonator

Q_static = Qmax_static;
static_resonator_width = f_0 / Q_static;
Q_dynamic = Qmax_dynamic;
dynamic_resonator_width = f_0 / Q_dynamic;

f = f_0 * ones(size(t));

static_signal = resonator(signal, f_0, f, Q_static, f_sampling);
static_noise  = resonator(noise,  f_0, f, Q_static, f_sampling);
static_gain = bandpower(static_signal) / bandpower(signal);
static_SNR = bandpower(static_signal) / bandpower(static_noise)

%display_spectra(signal, static_signal, f_sampling);
%hold on
%{
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
%}
%uiwait
f = f_iwave;
dynamic_signal = resonator(signal, f_0, f, Q_dynamic, f_sampling);
dynamic_noise  = resonator(noise,  f_0, f, Q_dynamic, f_sampling);
%dynamic_gain = bandpower(dynamic_signal) / bandpower(signal);
dynamic_SNR = bandpower(dynamic_signal) / bandpower(dynamic_noise)
%display_spectra(static_noise, static_signal, f_sampling);

improvement = dynamic_SNR / static_SNR
end
