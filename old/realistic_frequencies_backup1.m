f_sampling = 1000e6; % Hz
t_end = 0.002; % secs

t_step = 1/f_sampling;
t = 0 : t_step : t_end;

num_timesteps = length(t);

%% uncorrelated gaussian points
f_uncorrelated = randn(1, num_timesteps);

%% correlate the points
correlation_time = 0.0004; % secs
assert(correlation_time < t_end / 3) % this not being true caused a long problem once
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

figure
subplot(2, 1, 1)
[acf, k] = autocorr(f_uncorrelated, NumLags=num_timesteps-1);
plot(k/f_sampling, acf)
subplot(2, 1, 2)
[acf, k] = autocorr(f_correlated, NumLags=num_timesteps-1);
plot(k/f_sampling, acf)
yline(0)
return


% renormalise f_correlated to have std.dev = 1
f_renormalised = f_correlated / std_dev_correlated;
% confirm this worked
std_dev_renormalised = std(f_renormalised);

f_0 = 1e6; % Hz, default axion frequency
axion_linewidth = f_0 * 1e-6; %1000 * f_0 * (240/3e5)^2; % Hz, spread of axion signal

% array in time of correlated frequencies with given linewidth
f = f_0 + axion_linewidth * f_renormalised;

% calculate a wave in time with the given frequency
x = calc_wave_with_frequency(t, t_step, f);

small_extract = 1e5 : 1.005e5;

% calculate avg df/dt (for calculating tau_iwave)
dfdsample = f(1:end-1) - f(2:end);
dfdt = dfdsample * f_sampling;
dfdt_square = dfdt .* dfdt;
dfdt_square_avg = mean(dfdt_square);

tau_iwave = (288 * pi^4 * f_0^2 * dfdt_square_avg )^(-1/6); % tau_opt from the iwave paper
%tau_iwave = 5e-4;

f_iwave = calc_iwave(f_sampling, f_0, t, x, tau_iwave);
% plot iwave before and after
%{
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

% now we feed our wave into a resonator tuned at (a) the default axion
% frequency and (b) the moving frequency from the iwave calculation

%y = zeros(size(x));
%w = 1;
%    % start at 2nd element
%for i = 2 : length(x) % should it be Delta_0(i) or (i-1)?
%    y(i) = exp(-w + 1i * Delta_0(i))*y(i-1) + (1 - exp(-w))*x(i-1);
%end

x = randn(1, num_timesteps);
%{
carrier_frequency = 8.5*f_0;

osc      =  cos(2*pi*carrier_frequency*t);
osc_quad = -sin(2*pi*carrier_frequency*t);

I = x .* osc;
Q = x .* osc_quad;

x_q = I + 1i * Q;

y_q = zeros(size(x_q));
w = 0.1;
Delta_0 = 1e3; %*ones(size(x_q));
    % start at 2nd element
for i = 2 : length(x_q) % should it be Delta_0(i) or (i-1)?
    y_q(i) = exp(-w + 1i * Delta_0)*y_q(i-1) + (1 - exp(-w))*x_q(i-1);
end

I = real(y_q) .* osc;
Q = imag(y_q) .* osc_quad;

y = I + Q;

display_spectra(x, y, f_sampling)
title('straight code')
uiwait
%}

%y = resonator_ed(x, t, 1e9, 5e3, f_sampling);
%display_spectra(x, y, f_sampling)
%uiwait

f_0 = 1e9;
Q = 5e3;

y = resonator_ed(x, t, f_0, Q, f_sampling);
f = display_spectra(x, y, f_sampling);
hold on

%s = 1i * 2 * pi * f;
%omega_res = 2*pi* 2.5e6;
%width = 1e7;
%Hs = 0.4 * omega_res^2 ./ ( s.^2 + width*s + omega_res^2 );
%splot(f, abs(Hs), 'DisplayName','H(s) mine')

omega_0 = 2*pi * f_0;
Delta_0 = omega_0 / f_sampling;
w = Delta_0 / (2 * Q);
Delta = 2*pi*f/f_sampling;
A_Delta = (1 - exp(-w)) ./ sqrt(1 - 2*exp(-w)*cos(Delta - Delta_0) + exp(-2*w)); % from Ed paper
plot(f, A_Delta, 'DisplayName','A(Delta)')

%omega_0 = 2*pi * 2e6;
%Gamma = omega_0 / Q;
%Hs = Gamma * s ./ ( s.^2 + Gamma*s + omega_0^2 );
%plot(f, abs(Hs), 'DisplayName','H(s) Ed')

uiwait

%{
ps = 10 .^ (3:8)
fs = zeros(size(p));

for i = 1:6
    w = 0.1;
    Delta_0 = ps(pow)

    y = resonator_ed(x, t, Delta_0, w);

    fs(i) = get_ratio_f_0(x, y, f_sampling);
    display_spectra(x, f, f_sampling)
    uiwait
end

scatter(ps, fs)
%}
%{
plot(t(small_extract), x(small_extract))
hold on
plot(t(small_extract), y(small_extract))
legend('x', 'y')
uiwait
%}

%{
%x = randn(1, num_timesteps);

Q_resonator = 1e3;

%y_tracking = resonator(x, f_iwave - 5e5, Q_resonator, f_sampling, num_timesteps);
%figure(1);
%display_spectra(x, y_tracking, f_sampling)

%f_sweep = linspace(f_0/2, f_0*2, num_timesteps);
%x = calc_wave_with_frequency(t, t_step, f_sweep);
f_static_val = f_0*10;
f_static = f_static_val * ones(1, num_timesteps);
y_static = resonator(x, f_static, Q_resonator, f_sampling, num_timesteps);
%figure(2);
display_spectra(x, y_static, f_sampling)

bandpass_centre = f_0;
bandpass_width = 1e5;

band = [bandpass_centre - bandpass_width/2, bandpass_centre + bandpass_width/2];
%y_tracking = bandpass(y_tracking, band, f_sampling);
y_static   = bandpass(y_static,   band, f_sampling);

%power_tracking = rms(y_tracking)
power_static   = rms(y_static)
% maybe should use bandpower() instead of rms()

% plot theoretical response
hold on
f_0_theoretical = f_static_val;
Q_theoretical = Q_resonator;
f = 1e4 : 1e4 : 1e8;
omega = 2*pi*f;
s = 1i * omega;
Gamma = 2*pi*f_0_theoretical / Q_theoretical; % 1e3 is Q
omega_0 = 2*pi*f_0_theoretical;
H = Gamma.*s ./ ( s.^2 + Gamma.*s + omega_0^2);
loglog(f, abs(H), 'DisplayName','H(s)')
%}

%{
Gc = dBtoRatio(130);
Gw = dBtoRatio(0); % is this right?
A = dBtoRatio(-130);
assert(A*Gc*Gw <= 1, "see eq 4 of Daw 2019")
FWHM = 1e6 * 1e-6; % copied from axion_linewidth
Gamma = 2*pi * FWHM;
omega_0 = 2*pi * 1e6;
extra = (1+A*Gc*Gw) / Gc;
%y = filter([0, extra * Gc * Tau], [1, (1 + A*Gc*Gw)*Tau, omega_0^2],x);
y_tracking = filter([0, omega_0], [1, omega_0], x);
%}
%{
[amp_bef, f] = nsd_pwelch(x, 0.002, f_sampling);
[amp_aft, f] = nsd_pwelch(y, 0.002, f_sampling);
subplot(3,1,1)
loglog(f, amp_bef)
subplot(3,1,2)
loglog(f, amp_aft)
subplot(3,1,3)
ratio = amp_aft./amp_bef;
loglog(f, ratio)
uiwait
%}

%{
f = 1e1 : 1e4 : 1e9;
omega = 2*pi*f;
s = 1i * omega;
H = ( extra * Gc*Tau.*s ) ./ ( s.^2 + (1 + A*Gc*Gw)*Tau.*s + omega_0^2);
loglog(f, abs(H))
%hold on
%uiwait
%}

function y = resonator_ed(x, t, f_0, Q, f_sampling)
    % this all from Daw 2019 appendix
    % atm this is actually a resonant low pass filter

    omega_0 = 2*pi * f_0;
    Delta_0 = omega_0 / f_sampling;
    w = Delta_0 / (2 * Q);  

    carrier_frequency = 2.5e6;
    osc      =  cos(2*pi*carrier_frequency*t);
    osc_quad = -sin(2*pi*carrier_frequency*t);

    I = x .* osc;
    Q = x .* osc_quad;
    
    x_q = I + 1i * Q;
    y_q = zeros(size(x_q));

        % start at 2nd element
    for i = 2 : length(x_q) % should it be Delta_0(i) or (i-1)?
        y_q(i) = exp(-w + 1i * Delta_0)*y_q(i-1) + (1 - exp(-w))*x_q(i);
    end
    
    I = real(y_q) .* osc;
    Q = imag(y_q) .* osc_quad;
    
    y = I + Q;
end

function y = resonator(x, f, Q, f_sampling, num_timesteps)
    % method using matlab builtins, doesn't support variable centre freq
    %[b, a] = iirpeak(2*f(1)/f_sampling, 1/Q);
    %y = filter(b, a, x);
    
    y = zeros(1, num_timesteps);
        % start at 3rd
    for i = 3 : num_timesteps
        %f_resonator = 2e5;
        f_resonator = f(i);
        Q_resonator = Q;

        [b, a] = my_iirpeak(f_resonator, Q_resonator, f_sampling);
        
        %wo = f_resonator / (f_sampling/2);
        %bw = wo / Q_resonator;
        %[b, a] = iirpeak(wo, bw);
        y(i) = b(1)*x(i) + b(2)*x(i-1) + b(3)*x(i-2) - a(2)*y(i-1) - a(3)*y(i-2);
    end
end

function f = display_spectra(x, y, f_sampling)
    [amp_bef, f] = nsd_pwelch(x, 0.0005, f_sampling);
    [amp_aft, f] = nsd_pwelch(y, 0.0005, f_sampling);
    subplot(2,1,1)
    loglog(f, amp_bef)
    hold on 
    loglog(f, amp_aft)
    legend('before', 'after')
    %lim_x = xlim();
    %lim_y = [1e-10, 1];
    %ylim(lim_y)
    %xlim(lim_x)
    %subplot(3,1,2)
    %ylim(lim_y)
    %xlim(lim_x)
    subplot(2,1,2)
    ratio = amp_aft./amp_bef;
    loglog(f, ratio)
    legend('transfer')
    [ratio_max, index_max] = max(ratio);
    f_0_ratio = f(index_max)
    %ylim(lim_y)
    %xlim(lim_x)
end

function ratio = dBtoRatio(dB)
    ratio = 10^(dB / 20); % or dB/10?
end

function x = calc_wave_with_frequency(t, t_step, f)
    x = zeros(size(t));
    f_integral = 0;

    for i = 1 : length(t)
        x(i) = sin(2*pi*f_integral);
        f_integral = f_integral + f(i)*t_step;
    end
end

function freq = calc_iwave(f_sampling, f_guess, time, wave, tau)
    % initialise iwave pll
    start_time = 0.0; % secs
    fstate = iwavePllConstruct(f_sampling, tau, f_guess, start_time);
    
    % run iwave pll
    [dout, qout, amp, error, freq, fstate] = iwavePllRun(wave, fstate);
end

function [b, a] = my_iirpeak(f_0, Q, f_sampling)
    % lifted straight from scipy source code: https://github.com/scipy/scipy/blob/v1.8.0/scipy/signal/_filter_design.py#L4946-L5023
    w0 = f_0 / (f_sampling / 2); % normalise to nyquist frequency

    assert(w0 > 0)
    assert(w0 < 1)

    w0 = w0*pi;

    % Compute beta: formula 11.3.19 (p.579) from reference [1]
    beta = tan(w0/(2*Q));

    % Compute gain: formula 11.3.6 (p.575) from reference [1]
    gain = 1 / (1 + beta);

    % Compute numerator b and denominator a
    % formulas 11.3.7 (p.575) and 11.3.21 (p.579)
    % from reference [1]
    b = (1 - gain) .* [1, 0, -1];
    a = [1.0, -2*gain*cos(w0), 2*gain - 1];
end