function f_axion = generate_axion_frequency(t, f_sampling, correlation_time, f_0, axion_linewidth)

%% uncorrelated gaussian points
f_uncorrelated = randn(size(t));

%% correlate the points
assert(correlation_time < t(end) / 2) % this not being true caused a long problem once
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

%% turn into an array of frequencies
% renormalise f_correlated to have std.dev = 1
f_renormalised = f_correlated / std_dev_correlated;
% confirm this worked
std_dev_renormalised = std(f_renormalised);
% now do the axion stuff, commented out while we work just on noise

% array in time of correlated frequencies with given linewidth
f_axion = f_0 + axion_linewidth * f_renormalised;

%f_axion = linspace(4e7, 6e7, length(t));

end