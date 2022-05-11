m1 = -1.8;
m2 = 0.96;
% Hz = (z + 1).^2 / (z.^2 + m1*z + m2);
t_step = 1;
omega = logspace(-2, 0, 1e4);
z = exp(t_step * 1i * omega);
Hw = (z + 1).^2 ./ (z.^2 + m1*z + m2);
s = 1i * omega
Hs = 25 ./ (s.^2 + 0.1*s + 0.15);
loglog(omega, abs(Hw));
hold on
loglog(omega, abs(Hs));
legend('Hw', 'Hs')