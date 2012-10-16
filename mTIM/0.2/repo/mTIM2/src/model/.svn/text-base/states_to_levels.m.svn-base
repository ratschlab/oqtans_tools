function level = states_to_levels(state_seq, PAR)

% level = states_to_levels(state_seq, PAR)
%
% Parses the predicted expression level from a state sequence. Non-exonic
% states will have level 0.
%
% state_seq -- sequence of states of lenth n
% PAR -- a struct of parameters specified in setup_training.m and
%   train_hmsvm.m
% returns a sequence of predicted expression levels of length n (which is
%   0 for al probes not predicted to be exonic, and between 1 and
%   PAR.num_levels else)
%
% written by Georg Zeller, MPI Tuebingen, Germany 2008-2009

% get label set and state names
LABELS = get_label_set_mTIM();
STATES = get_state_set_mTIM(PAR);

% get state ids of exon states (soreted by level)
fn = fieldnames(STATES);
exw_states = [strmatch('EFW', fn), strmatch('EIW', fn), strmatch('ELW', fn)];
exc_states = [strmatch('EFC', fn), strmatch('EIC', fn), strmatch('ELC', fn)];
assert(size(exw_states,1) == PAR.num_levels);
assert(size(exc_states,1) == PAR.num_levels);
num_ex_state_sets = size(exw_states,2);
assert(num_ex_state_sets == size(exc_states,2));

% convert state sequence into level sequence
level = zeros(size(state_seq)); % nonexonic states will have level 0
for s=1:num_ex_state_sets,
  for i=1:PAR.num_levels,
    idx = find(state_seq == exw_states(i,s) | state_seq == exc_states(i,s));
    level(idx) = i;
  end
end
% eof