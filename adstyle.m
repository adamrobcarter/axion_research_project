function adstyle(width,height)
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
figure

set(gcf,'DefaultAxesFontName','Helvetica');
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontWeight','normal');

set(gcf,'PaperType','a4');
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');

set(gcf,'units','centimeters','position',[16, 8, width, height])

box on
axis padded
