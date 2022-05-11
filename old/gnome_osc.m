f_sampling = 1024; % Hz
t_end = 1; % secs

t_step = 1/f_sampling;
t = 0 : t_step : t_end;

num_timesteps = length(t);

%% uncorrelated gaussian points
f_uncorrelated = randn(1, num_timesteps);

%% correlate the points
correlation_time = 0.1; % secs
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
histogram(f)
subplot(2, 1, 2)
histogram(f_correlated)
%}

% the moving mean averaging clearly reduces the std. deviation of the
% points, but by how much? must be a theoretical answer, or could just
% calculate from the data. By inspection, it appears to be reduced by a 
% factor of sqrt(boxcar_length)
std_dev_gaussian = std(f_uncorrelated);
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
%}

% renormalise f_correlated to have std.dev = 1
f_renormalised = f_correlated / std_dev_correlated;
% confirm this worked
std_dev_renormalised = std(f_renormalised);


f_0 = 100; % Hz, default axion frequency
axion_linewidth = f_0 / 10; % spread of axion signal

% array in time of correlated frequencies with given linewidth
f = f_0 + axion_linewidth * f_renormalised;

wave = calc_wave_with_frequency(t, t_step, f);

% calculate avg df/dt (for calculating tau_iwave)
dfdsample = f(1:end-1) - f(2:end);
dfdt = dfdsample * f_sampling;
dfdt2 = dfdt .* dfdt;
dfdt2_avg = mean(dfdt2);

tau_iwave = (288 * pi^4 * f_0^2 * dfdt2_avg )^(-1/6)

f_iwave = calc_iwave(f_sampling, f_0, t, wave, tau_iwave);

figure
plot(t, f)
hold on
plot(t, f_iwave)
legend('before iwave', 'after iwave')
xlabel('time (s)')
ylabel('frequenzy (Hz)')

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

    % plot results
    %{
    figure;
    subplot(2,2,1);
    plot(time,amp,'-');
    xlabel('time (s)');
    ylabel('amplitude');
    subplot(2,2,2);
    plot(time,error,'-');
    xlabel('time (s)');
    ylabel('error');
    subplot(2,2,3);
    plot(time,freq,'-');
    xlabel('time (s)');
    ylabel('frequency (Hz)');
    subplot(2,2,4);
    plot(time,dout,'-',time,qout,'-');
    xlabel('time (s)');
    ylabel('D and Q out');
    
    uiwait
    %}
end
%}