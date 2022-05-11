function fig=edstyle(leftmargin,bottommargin,width,height)
% EDSTYLE - make a matlab figure which exports to a 
% publication-ready file.
%
% useage: <figure handle> = edstyle(<left margin (/a4)>, ...
%                                   <bottom margin (/a4)>, ...
%                                   <width> (/a4), ...
%                                   <height> (/a4))
%
% Ed Daw, 13th January 2014
%     rev 3rd Feburary 2014
 
% make plot size defaults
set(0,'DefaultAxesFontName','Helvetica');
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultAxesFontWeight','normal');
fig=figure;
framesize=[leftmargin bottommargin width height];
set(fig,'PaperType','a4');
set(fig,'PaperPositionMode','manual');
set(fig,'PaperUnits','centimeters');
set(fig,'PaperPosition',framesize);

