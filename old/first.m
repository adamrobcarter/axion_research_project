f_sampling = 1e3;
t_end = 200;

t = 0 : 1/f_sampling : t_end;

f = 100;

omega = 2*pi*f;

x = cos(omega*t);

WinTime = 1;
    
nsd_pwelch(x, 1, f_sampling);