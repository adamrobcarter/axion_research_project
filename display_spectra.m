function [ratio, f] = display_spectra(x, y, f_sampling, range)
arguments
    x
    y
    f_sampling
    range = [0, f_sampling/2]
end
% show two spectra and their transfer function

[amp_bef, f] = nsd_pwelch(x, 0.001, f_sampling);
[amp_aft, f] = nsd_pwelch(y, 0.001, f_sampling);

fig = figure;
%set(fig,'Units','normalized','Position',[0 0 1 .7]); 

subplot(2,1,1)
loglog(f, amp_bef)
hold on 
loglog(f, amp_aft)
legend('before', 'after', 'Location', 'southwest')
%xlim([1e7, f_sampling/2])
title('simulated power spectrum before/after resonant filter')
xlim(range)

subplot(2,1,2)
ratio = amp_aft./amp_bef;
loglog(f, ratio)
legend('simulation', 'Location', 'southwest')
title('filter transfer function')
xlim(range)

[pks,locs,w,p] = findpeaks(ratio, WidthReference="halfheight", SortStr='descend');
maxindex = 1;
peak_freq = f(locs(maxindex));
peak_w = w(maxindex);
%peak_Q = peak_freq / peak_w

end