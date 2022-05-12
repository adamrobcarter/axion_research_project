function improvement = get_improvement(signal, noise, f_0, f_sampling, f_axion, f_iwave, t)

avg_dist_dynamic = rms(f_iwave - f_axion);
Qmax_dynamic = 2 * f_0 / avg_dist_dynamic;
avg_dist_static = rms(f_0 * ones(size(t)) - f_axion);
Qmax_static = 2 * f_0 / avg_dist_static;

Q_static = Qmax_static;
%static_resonator_width = f_0 / Q_static;
Q_dynamic = Qmax_dynamic;
%dynamic_resonator_width = f_0 / Q_dynamic;

f = f_0 * ones(size(t));

static_signal = resonator(signal, f_0, f, Q_static, f_sampling);
static_noise  = resonator(noise,  f_0, f, Q_static, f_sampling);
%static_gain = bandpower(static_signal) / bandpower(signal);
static_SNR = bandpower(static_signal) / bandpower(static_noise);

f = f_iwave;
dynamic_signal = resonator(signal, f_0, f, Q_dynamic, f_sampling);
dynamic_noise  = resonator(noise,  f_0, f, Q_dynamic, f_sampling);
%dynamic_gain = bandpower(dynamic_signal) / bandpower(signal);
dynamic_SNR = bandpower(dynamic_signal) / bandpower(dynamic_noise);
%display_spectra(static_noise, static_signal, f_sampling);

improvement = dynamic_SNR / static_SNR;

end

