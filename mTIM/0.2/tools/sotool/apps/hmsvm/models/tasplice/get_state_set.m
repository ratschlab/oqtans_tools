function STATES = get_state_set(PAR)
% STATES = get_state_set(PAR)
% returns the set of states of the graphical model
%  as well as the number of discrete expression levels
%
% written by Georg Zeller, MPI Tuebingen, Germany

STATES = [];
STATES.ige         =  1; % intergenic
STATES.ige_ss      =  2; % splice site state between intergenic probes
STATES.trans_start =  3; % transcript start
STATES.trans_end   =  4; % transcript end

STATES.exo         =  5; % exonic
STATES.exo_ss      =  6; % splice site state between adjacent probes of
                         % the same exon
STATES.don         =  7; % splice site state between exo -> ino
STATES.acc         =  8; % splice site state between ino -> exo
STATES.ino         =  9; % intronic
STATES.ino_ss      = 10; % splice site states between adjacent probes of
                         % the same intron
