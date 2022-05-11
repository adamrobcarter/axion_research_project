addpath 'C:\Users\anonf\OneDrive\uni\phys project\matlab\iwave\development\matlab'

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
axion_power = 1e-21; % W, from eq 2.28 of Daw thesis

% noise powers 1e-16 2e-16

%% axion signal params
f_0 = axion_rest_energy / h_eV; % Hz, default axion frequency
axion_linewidth = f_0 * 0.5 * (axion_velocity / c)^2;
axion_Q = f_0 / axion_linewidth;
de_broglie_wavelength = 2*pi * hbar_times_c / (axion_rest_energy * axion_velocity/c);
t_coherence = de_broglie_wavelength / axion_velocity;

%% generate axion frequency
f_axion = generate_axion_frequency(t, f_sampling, t_coherence, f_0, axion_linewidth);

%% generate signal and noise
signal = calc_wave_with_frequency(t, f_sampling, f_axion, axion_power);
%noise_power = calculate_noise_power(f_sampling, T);
noise = randn(size(t)); % this noise has total power 1, now scale it

noise_powers = 1 : 5;
noise_powers = 1e-16 * noise_powers
reps = 1;

tau_res = zeros([length(noise_powers), reps]);

for i = 1 : length(noise_powers)
    for rep = 1 : reps
        tau_res(i, rep) = get_tau(noise_powers(i), f_sampling, t, signal, noise, f_0, f_axion)
    end
end

function tau_iwave = get_tau(noise_power, f_sampling, t, signal, noise_in, f_0, f_axion)

    noise = noise_in * sqrt(noise_power);
    
    % atm we're running iwave on the signal not the signal + noise
    iwave_input = signal + noise;
    
    band_width = 150;
    band = [f_0-band_width/2, f_0+band_width/2];
    iwave_input = bandpass(iwave_input, band, f_sampling);
    
    % front to chop off
    start_secs = t(end)/4;
    start = start_secs * f_sampling;
    assert(start_secs < t(end)/3)
    
    chop = @(arr) arr(start:end);
    to_minimise = @(tau) sum(power(chop(f_axion - run_iwave(iwave_input, f_sampling, f_0, tau)), 2));
    
    %options = optimoptions(@fmincon, 'MaxFunctionEvaluations', 12, 'Display','iter');
    options = optimoptions(@fmincon, 'MaxIterations', 4, 'Display','iter');
    [x,fval,exitflag,output] = fmincon(to_minimise, 5e-4, [], [], [], [], 1/f_sampling, t(end), [], options);
    tau_iwave = x;
end

%f_iwave = run_iwave(iwave_input, f_sampling, f_0, tau_iwave);

% plot iwave before and after
%{
figure
plot(t, f_axion)
%lim = ylim();
hold on
plot(t, f_iwave)
%ylim(lim)
legend('before iwave', 'after iwave')
xlabel('time (s)')
ylabel('frequenzy (Hz)')
xline(t_coherence, '-', 't_c')
%xline(tau_iwave, '-', 'tau')
xline(start_secs, '-', 'chop')
%}