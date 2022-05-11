function x = calc_wave_with_frequency(t, f_sampling, f)
% calculate a wave in time with the given frequency in time

t_step = 1 / f_sampling;

x = zeros(size(t));
f_integral = 0;

for i = 1 : length(t)
    x(i) = sin(2*pi*f_integral);
    f_integral = f_integral + f(i)*t_step;
end

end