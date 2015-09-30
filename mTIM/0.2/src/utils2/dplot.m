function [ax1, ax2]=dplot(x1, y1, lprop1, axprop1, axprop2) ;
% [ax1, ax2]=dplot(x1, y1, prop1, prop2) ;

%clf ;
hl1=line(x1, y1, lprop1) ;
ax1=gca ;
set(ax1, axprop1) ;

ax2=axes('Position', get(ax1, 'Position'), 'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', 'Color', 'none', 'YLim', get(ax1, 'YLim'), ...
    axprop2) ;
