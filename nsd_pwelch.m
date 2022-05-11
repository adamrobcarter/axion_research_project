function [V,f] = nsd_pwelch(v,WinTime,Rate)
% NSD_PWELCH - generate amplitude spectral density of a data stream
%
% useage: [V,f]=nsd_pwelch(x,WinTime,Rate), where
% f is an array of frequencies
% V is an array of amplitude spectral densities at those frequencies
% x is the input data array
% WinTime is the time duration of each short Fourier transform used
%     to form the Welch estimate. Note that the bandwidth of each bin
%     of the FFT is the reciprocal of this window time. Longer times
%     lead to higher resolution averages, so there are more bins, but
%     there are less FFTs to average; consequently the bin-to-bin noise
%     fluctuations are larger. It's a trade-off.
% Rate is the sampling rate in Hz
%
% Elena Massera and Ed Daw (e.daw@sheffield.ac.uk), 2017/04/04, based
% on the nsd() function from the MIT group, Shuorov Chatterji, 1998.

WinLen = 2*round((Rate * WinTime)/2);
FFTLen = 2 ^ ceil( log(WinLen) / log(2) );
Fs=Rate;
N=FFTLen;
Data=v;

[PwrSpc,f]= pwelch(Data, hanning(N),0, N, Fs);
V=sqrt(PwrSpc);

% make plot if the output arguments are not given in the call
if ~nargout,
  loglog(f,V);    
    xlabel('frequency [Hz]');
    ylabel('amplitude noise (amplitude/\surdHz)');
    grid;
end
