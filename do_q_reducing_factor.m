function signal_mod = do_q_reducing_factor(signal, f_axion, f_0, Q)
% f_0 can be single number or array of size(signal)=size(f_axion)

factor = 1 ./ ( 1 + 4 * Q^2 * (f_axion - f_0).^2 ./ f_0.^2 );

signal_mod = signal .* factor;

end