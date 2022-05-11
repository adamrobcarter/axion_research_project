x = randn(1, 1e8);
[a, f] = nsd_pwelch(x, 0.001, 1e8);
loglog(f, a);