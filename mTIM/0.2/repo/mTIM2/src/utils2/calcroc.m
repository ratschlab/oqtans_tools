function [false_alarms,hits] = calcroc(output,LTE)

assert(all(size(output)==size(LTE))) ;

[ld,idx] = sort(output); 

i=0 ;

hits=1-cumsum(LTE(idx)>0)/sum(LTE > 0) ;
false_alarms=1-cumsum(LTE(idx)<0)/sum(LTE < 0) ;

% Da die Kurve nie alle Negativen als positiven zaehlt, hat
% Falsealarms nie 
% den Wert 1. Cumsum zaehlt das erste Element als negatives, diese
% fehlt bei
%  den False_alarms.
% Daher wird vorne an den Vektor noch eine 1 gehaengt

hits = [1 hits ];
false_alarms = [1 false_alarms];

if nargout<2,
  false_alarms = [false_alarms;hits];
end ;
