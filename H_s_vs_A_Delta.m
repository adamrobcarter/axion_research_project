f_sampling = 1e8;
f = logspace(1, log10(f_sampling/2), 1e4);

f_0 = 1e5;
Q = 5e3;

omega_0 = 2*pi * f_0;
Delta_0 = omega_0 / f_sampling;
w = Delta_0 / (2 * Q);
Delta = 2*pi*f/f_sampling;
A_Delta = (1 - exp(-w)) ./ sqrt(1 - 2*exp(-w)*cos(Delta - Delta_0) + exp(-2*w)); % from Ed paper
loglog(f, A_Delta, 'DisplayName','A(Delta)')

hold on

Gamma = omega_0 / Q;
s = 1i * 2 * pi * f;
Hs = Gamma * s ./ ( s.^2 + Gamma*s + omega_0^2 );
loglog(f, abs(Hs), 'DisplayName','H(s)')

legend()