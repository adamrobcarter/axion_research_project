f_sampling = 1000e6; % Hz
t_end = 0.002; % secs

t_step = 1/f_sampling;
t = 0 : t_step : t_end;

num_timesteps = length(t);

%% uncorrelated gaussian points
f_uncorrelated = randn(1, num_timesteps);

%% correlate the points
correlation_time = 0.00001; % secs
assert(correlation_time < t_end / 10) % this not being true caused a long problem once
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

x = randn(1, num_timesteps);

f_0 = 1e6;
Q = 5e3;

y = resonator(x, f_0, Q, f_sampling);
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
Hs = Gamma * s ./ ( s.^2 + Gamma*s + omega_0^2 );
plot(f, abs(Hs), 'DisplayName','|H(s)|')

uiwait

function y = resonator(x, f_0, Q, f_sampling)
    % this all from Daw 2019 appendix
    % atm this is actually a resonant low pass filter

    omega_0 = 2*pi * f_0;
    Delta_0 = omega_0 / f_sampling;
    w = Delta_0 / (2 * Q);  
    
    I = x;
    quarter_period_secs = 1/f_0 * 1/4;
    quarter_period_samples = quarter_period_secs * f_sampling;
    Q = delay(x, quarter_period_samples);
    
    x_q = I + 1i * Q;
    y_q = zeros(size(x_q));

    for i = 2 : length(x_q) % should it be Delta_0(i) or (i-1)?
        y_q(i) = exp(-w + 1i * Delta_0)*y_q(i-1) + (1 - exp(-w))*x_q(i);
    end

    I = delay(real(y_q), quarter_period_samples);
    Q = imag(y_q);
    
    y = I + Q;
end

% delay a signal by the specified number of samples.
% pads the start with zeroes
function x_delayed = delay(x, delay_samples)
    x_delayed = zeros(size(x));
    x_delayed(delay_samples+1:end) = x(1:end-delay_samples);
end

function f = display_spectra(x, y, f_sampling)
    [amp_bef, f] = nsd_pwelch(x, 0.0005, f_sampling);
    [amp_aft, f] = nsd_pwelch(y, 0.0005, f_sampling);

    subplot(2,1,1)
    loglog(f, amp_bef)
    hold on 
    loglog(f, amp_aft)
    legend('before', 'after')

    subplot(2,1,2)
    ratio = amp_aft./amp_bef;
    loglog(f, ratio)
    legend('simulation')
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