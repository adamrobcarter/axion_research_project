function [pdout, pqout, paout, peout, pfout, pfstateout] = iwave_my_run(fs_in, tau_in, fline_in, starttime, INPUT_DATA)

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

% translated from iwave/development/src/iwave_pll_run.c

datalength = length(INPUT_DATA);
pindata = INPUT_DATA
pdout = zeros(1, datalength);
pqout = zeros(1, datalength);
paout = zeros(1, datalength);
peout = zeros(1, datalength);
pfout = zeros(1, datalength);
pfstateout = zeros(1, datalength);

%foundnan=0;
%double yr,yi,yrm,yim,yr2,yi2,yrm2,yim2;
%double x2r,x2i,asq,esig,delta,deltausb,srate;
%double w,absq;
%double deltaprefold;


% load the previous output elements
yrm=ps.ymemreal;
yim=ps.ymemimag;
yrm2=ps.ymemrealusb;
yim2=ps.ymemimagusb;

% run the filter over the input data
for dcount = 1 : datalength
    % trap NANs
    %if isnan(pindata(dcount))!=0
    %    pindata(dcount)=0;
    %    foundnan=1;
    %    % show error here 
    %end
    % run 1f iwave
    yr=ps.adiag*yrm - ps.aoffdiag*yim + ps.b*pindata(dcount);
    yi=ps.aoffdiag*yrm + ps.adiag*yim;
    pdout(dcount) = ps.r11*yr + ps.r12*yi;
    pqout(dcount) = ps.r12*yr + ps.r22*yi;
    yrm=yr;
    yim=yi;
    % calculate phase discriminant
    asq=pdout(dcount)*pdout(dcount)+pqout(dcount)*pqout(dcount);  
    x2r=pindata(dcount)*pdout(dcount)+pqout(dcount)*pqout(dcount)-asq;
    x2i=(pindata(dcount)-pdout(dcount))*pqout(dcount);
    % run 2f iwave on the phase discriminant */
    yr2=ps.adiagusb*yrm2 - ps.aoffdiagusb*yim2 + ps.busb*x2r;
    yi2=ps.aoffdiagusb*yrm2 + ps.adiagusb*yim2 + ps.busb*x2i;
    yrm2=yr2;
    yim2=yi2;
    % derive the error signal
    if asq>0
      esig=2*(x2i-yi2)/asq;
    else
      esig=2*(x2i-yi2);
    end
    paout(dcount)=sqrt(asq);
    peout(dcount)=esig;
    % close loop by modifying delta
    deltaprefold = ps.delta - ps.gain*esig;
    % run single frequency foldback to stay within nyquist band
    delta=iwave_firstbz(deltaprefold);
    % write frequency to data output
    pfout(dcount) = delta*ps.srate / (2*pi);
    % recalculate the elements of the state data that are affected
    % by the change in delta due to the closed feedback loop
    w = ps.w;
    absq = (1-exp(-w))*(1-exp(-w))/(1-2*exp(-w)*cos(2*delta)+exp(-2*w));
    ps.adiag       = exp(-w)*cos(delta);
    ps.aoffdiag    = exp(-w)*sin(delta);
    ps.r12         = ps.b * ps.adiag / ps.aoffdiag;
    ps.r22         = 4/(1-absq) - ps.r11;
    ps.ymemreal    = yrm;
    ps.ymemimag    = yim;
    % deltausb is the frequency where the upper sideband appears in
    % radians per sample, accounting for aliasing where appropriate
    srate=ps.srate;
    % account for possible out of band 2*delta. */
    deltausb=iwave_firstbz(2*delta);
    ps.adiagusb    = exp(-2*w)*cos(deltausb);
    ps.aoffdiagusb = exp(-2*w)*sin(deltausb);
    ps.ymemrealusb = yrm2;
    ps.ymemimagusb = yim2;
    ps.delta       = delta;
end

end

function out = iwave_firstbz(deltain)
    deltaped=deltain+pi;
    out = abs(deltaped-2*pi*floor(deltaped/(2*pi))-pi);
end