function label_seq = states_to_labels(state_seq, state_model)

% label_seq = states_to_labels(state_seq, state_model)
%
% Converts a state sequence into a label sequence.
%
% state_seq -- a sequence of states (see get_state_set_mTIM.m)
% state_model -- graphical model (see make_model_mTIM.m)
% returns a label_sequence reconstructed from the given sequence of
%   states
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2009

state_labeling = [state_model.label];
label_seq = convert_states2labels(state_seq, state_labeling);

% eof