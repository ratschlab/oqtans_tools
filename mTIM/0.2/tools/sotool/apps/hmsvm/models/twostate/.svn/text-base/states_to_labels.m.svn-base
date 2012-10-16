function label_seq = states_to_labels(state_seq, state_model)

% label_seq = states_to_labels(state_seq, state_model)
%
% Converts a state sequence into a label sequence.
%
% state_seq -- a sequence of states (see get_state_set.m)
% state_model -- graphical model (see make_model.m)
% returns a label_sequence reconstructed from the given sequence of
%   states
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2007-2008

label_seq = [state_model(state_seq).label];

% eof