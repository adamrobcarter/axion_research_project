
function y = resonator_2(x, f_0, band_width, f_sampling)
% pass the supplied wave x through a resonator
% given centre frequency f_0_centre
% resonant frequency (array same size as x) f_0,
% Q, and f_sampling

y = zeros(size(x));

% do the algorithm
for i = 3 : length(x)
    [b, a] = my_iirpeak2(f_0(i), band_width, f_sampling);
    y(i) = b(1)*x(i) + b(2)*x(i-1) + b(3)*x(i-2) - a(2)*y(i-1) - a(3)*y(i-2);
end

end