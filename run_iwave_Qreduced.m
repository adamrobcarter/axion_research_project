function [freq, amp] = run_iwave(signal, noise, f_axion, Q, f_sampling, f_guess, tau)
% run the iwave algorithm to extract the frequency from the given wave

% calculate avg df/dt (for calculating tau_iwave_opt)
%{
dfdsample = f_axion(1:end-1) - f_axion(2:end);
dfdt = dfdsample * f_sampling;
dfdt_avg = rms(dfdt);

tau_iwave_opt = (288 * pi^4 * f_0^2 * dfdt_avg^2 )^(-1/6); % tau_opt from the iwave paper
%}

factor = 1 / ( 1 + Q^2 * ())

wave = signal_mod + noise;

% initialise iwave pll
start_time = 0.0; % secs
fstate = iwavePllConstruct(f_sampling, tau, f_guess, start_time);

% run iwave pll
[dout, qout, amp, error, freq, fstate] = iwavePllRun(wave, fstate);

end