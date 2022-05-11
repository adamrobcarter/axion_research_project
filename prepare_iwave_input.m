function input = prepare_iwave_input(signal, f_0, bandwidth, f_sampling)

wo = f_0 / (f_sampling/2);  
bw = bandwidth / (f_sampling/2);
[b, a] = iirpeak(wo, bw);
input = filter(b, a, signal);
    
end

