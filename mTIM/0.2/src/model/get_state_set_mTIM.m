function STATES = get_state_set_mTIM(PAR)

% STATES = get_state_set_mTIM(PAR)
%
% Returns the set of states of the graphical model.
%
% PAR -- a struct of parameters specified in setup_training.m and
%   train_hmsvm.m
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2009

assert(isfield(PAR, 'num_levels'));
assert(PAR.num_levels > 0);

STATES = [];

cnt = 1;
STATES.IGE = cnt; % intergenic
cnt = cnt + 1;

% W and C denote the forward (Watson) and reverse (Crick) strand 
% of the genome, respectively
strands = ['W'; 'C'];

for s=1:length(strands),
  for i=1:PAR.num_levels,
    % transcription start site (first exon state)
    STATES = setfield(STATES, sprintf('EF%s_%02i', strands(s), i), cnt);
    cnt = cnt + 1;

    % internal exon states
    STATES = setfield(STATES, sprintf('EI%s_%02i', strands(s), i), cnt);
    cnt = cnt + 1;

    % transcription termination site (last exon state)
    STATES = setfield(STATES, sprintf('EL%s_%02i', strands(s), i), cnt);
    cnt = cnt + 1;
    
    % initial (first) intron states (Watson strand)
    STATES = setfield(STATES, sprintf('IF%s_%02i', strands(s), i), cnt);
    cnt = cnt + 1;
    
    % internal intron states
    STATES = setfield(STATES, sprintf('II%s_%02i', strands(s), i), cnt);
    cnt = cnt + 1;
    
    % terminal (last)intron states
    STATES = setfield(STATES, sprintf('IL%s_%02i', strands(s), i), cnt);
    cnt = cnt + 1;
  end
end

% eof