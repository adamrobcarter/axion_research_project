f_sampling = 250e6; % Hz
t_step = 1/f_sampling;

%% physical constants
c = 3e8; % m/s
h_eV = 4.136e-15; % eV
k_B = 1.38065e-23; % m2 kg s-2 K-1
h = 6.62607015e-34; % m2 kg / s
hbar_times_c = 2e-7; % eV m

%% parameters
axion_rest_energy = 2e-6; % eV/c2
axion_velocity = 240e3; % m/s
T = 10e-3; % K
axion_power = 1.0e-21; % W, from eq 2.28 of Daw thesis
axion_power_over_Q = axion_power / 7e4;

%% axion signal params
f_a_original = axion_rest_energy / h_eV; % Hz, default axion
axion_linewidth = f_a_original * 1e-6;
axion_Q = f_a_original / axion_linewidth;

%% heterodyne
heterodyning_ratio = 20;
heterodyning_shift = f_a_original * (heterodyning_ratio - 1)/heterodyning_ratio;
f_a = f_a_original - heterodyning_shift;

de_broglie_wavelength = 2*pi * hbar_times_c / (axion_rest_energy * axion_velocity/c);
t_coherence = de_broglie_wavelength / axion_velocity;

% KSVZ to DFSZ is divide power by 7.2