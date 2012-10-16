function label_seq = states_to_labels(state_seq, state_model)
% label_seq = states_to_labels(state_seq, state_model)
% converts a state sequence into label sequence

% written by Georg Zeller, MPI Tuebingen, Germany

label_seq = [state_model(state_seq).label];