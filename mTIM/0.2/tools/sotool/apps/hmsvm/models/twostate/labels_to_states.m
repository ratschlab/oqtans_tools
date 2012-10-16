function [state_seq, problems] = labels_to_states(label_seq, state_model, signal, PAR)

% state_seq = labels_to_states(label_seq, state_model, signal, PAR)
%
% Converts a label sequence into a state sequence.
%
% label_seq -- a sequence of labels (see get_label_set.m)
% state_model -- graphical model (see make_model.m)
% signal -- the feature matrix (sequence of observations) of size m x n
%   where m is equal to the number of features and n the combined length
%   of the training sequences
% PAR -- a struct of parameters specified in setup_hmsvm_training.m and
%   train_hmsvm.m
% returns the sequence of states corresponding to the given label
%   sequence
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008
problems = [];

LABELS = get_label_set();
STATES = get_state_set(PAR);

state_seq = repmat(STATES.negative, 1, length(label_seq));
% so far, label_seq and state_seq are identical
state_seq(label_seq==LABELS.positive) = STATES.positive;
state_seq(1) = STATES.start;
state_seq(end) = STATES.stop;