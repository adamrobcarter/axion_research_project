function input = prepare_iwave_input(signal, f_0, bandwidth, f_sampling)
% f_0 is the centre of the band we do the pre-filtering with

Q = f_0 / bandwidth;

input = resonator(signal, f_0, f_0*ones(size(signal)), Q, f_sampling);
    
end

