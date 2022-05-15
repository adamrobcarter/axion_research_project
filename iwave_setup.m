function [STATE_PS] = iwave_setup(fs_in, tau_in, fline_in, starttime)

% translated from iwave/development/src/iwave_pll_construct.c

gain=1/(2*tau_in*tau_in*fs_in*fs_in);
% initialise single-f static iwave
w=1/(fs_in*tau_in);
delta=2*pi*fline_in/fs_in;
absq = (1-exp(-w))*(1-exp(-w))/(1-2*exp(-w)*cos(2*delta)+exp(-2*w));
ps.adiag         = exp(-w)*cos(delta);  
ps.aoffdiag      = exp(-w)*sin(delta);
ps.b             = 1-exp(-w);
ps.r11           = 2 - ps.b;
ps.r12           = -ps.b * ps.adiag / ps.aoffdiag;
ps.r22           = 4/(1-absq) - ps.r11;
ps.ymemreal      = 0;
ps.ymemimag      = 0;
% deltausb is the frequency where the upper sideband appears in radians per sample
if(fline_in > (fs_in / 4))
  % account for aliasing when line frequency exceeds half nyquist
  deltausb=2*pi*fabs(2*fline_in-fs_in)/fs_in;
else
  % no need to account for aliasing for line frequency less than half nyquist
  deltausb=2*delta;
end
ps.adiagusb      = exp(-2*w)*cos(deltausb);
ps.aoffdiagusb   = exp(-2*w)*sin(deltausb);
ps.busb          = 1-exp(-2*w);
ps.ymemrealusb   = 0;
ps.ymemimagusb   = 0;
ps.delta         = delta;
ps.gain          = gain;
ps.w             = w;
ps.srate         = fs_in;

STATE_PS = ps;

end