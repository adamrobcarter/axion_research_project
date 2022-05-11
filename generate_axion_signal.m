function signal = generate_axion_signal(t, f_sampling, f_axion, axion_power)

%% calculate a wave in time with the given frequency
signal = calc_wave_with_frequency(t, f_sampling, f_axion);
signal_power = bandpower(signal); % = 0.5000 - must be a reason for this exactness
signal = signal / sqrt(signal_power) * sqrt(axion_power);
power_signal = bandpower(signal);

end