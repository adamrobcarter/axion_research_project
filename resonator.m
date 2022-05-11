
function y = resonator(x, f_0_centre, f_0, Q, f_sampling)
% pass the supplied wave x through a resonator
% given centre frequency f_0_centre
% resonant frequency (array same size as x) f_0,
% Q, and f_sampling

% this all from Daw 2019 appendix
% set up params
omega_0 = 2*pi * f_0;
Delta_0 = omega_0 / f_sampling;
w = Delta_0 / (2 * Q);  

% find the quadratures
I = x;
quarter_period_secs = 1/4 * 1/f_0_centre;
quarter_period_samples = quarter_period_secs * f_sampling;
Q = delay(x, quarter_period_samples);

x_q = I + 1i * Q;
y_q = zeros(size(x_q));

% do the algorithm
for i = 2 : length(x_q) % should it be Delta_0(i) or (i-1)?
    y_q(i) = exp(-w(i) + 1i * Delta_0(i))*y_q(i-1) + (1 - exp(-w(i)))*x_q(i);
end

% recombine the quadratures
I = delay(real(y_q), quarter_period_samples);
Q = imag(y_q);

y = I + Q;

end