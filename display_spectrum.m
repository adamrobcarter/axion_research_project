function [ratio, f] = display_spectrum(x, f_sampling, range)
arguments
    x
    y
    f_sampling
    range = [0, f_sampling/2]
end
% show two spectra and their transfer function

[amp_bef, f] = nsd_pwelch(x, 0.001, f_sampling);

fig = figure;
set(fig,'Units','normalized','Position',[0 0 1 .7]); 

%subplot(2,1,1)
loglog(f, amp_bef)
xlim(range)


end