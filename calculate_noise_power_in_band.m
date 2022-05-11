function noise_power = calculate_noise_power_in_band(T, f_bottom, f_top)

k_B = 1.38065e-23; % m2 kg s-2 K-1
h = 6.62607015e-34; % m2 kg / s

dPdf = @(f) 0.5*h*f + h*f/(exp(h*f/k_B/T)-1);
noise_power = integral(dPdf, f_bottom, f_top);

end