
function [state, y] = resonator_step(state, x, f_0, bandwidth, f_sampling)
% pass the supplied wave x through a resonator
% given centre frequency f_0_centre
% resonant frequency (array same size as x) f_0,
% Q, and f_sampling

[b, a] = my_iirpeak2(f_0, bandwidth, f_sampling);
y = b(1)*x + b(2)*state.x1 + b(3)*state.x2 - a(2)*state.y1 - a(3)*state.y2;

state.x2 = state.x1;
state.x1 = x;
state.y2 = state.y1;
state.y1 = y;

end