function noise_power = calculate_noise_power(f_sampling, T)

k_B = 1.38065e-23; % m2 kg s-2 K-1
h = 6.62607015e-34; % m2 kg / s

dPdf = @(f) 0.5*h*f + h*f/(exp(h*f/k_B/T)-1);
noise_power = integral(dPdf, 0, f_sampling/2);

end