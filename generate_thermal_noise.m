function noise = generate_thermal_noise(t, f_sampling, noise_power)

%% generate
noise = randn(size(t)); % this noise has total power 1, now scale it

noise = noise * sqrt(noise_power);
noise_power_meas = bandpower(noise);

% plot noise spectra
%f = f_0/2 : f_0*2
%hold on
%loglog(f, amp)
%loglog(f, dPdf(f))
%yline(power_noise)
%return

%SNR = power_signal / power_noise * sqrt( axion_linewidth * t_end);

end