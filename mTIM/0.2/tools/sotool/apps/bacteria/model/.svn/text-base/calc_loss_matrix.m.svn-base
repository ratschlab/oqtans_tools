function loss_matrix = calc_loss_matrix(true_state_seq, state_model, PAR)

% loss_matrix = calc_loss_matrix(true_state_seq, state_model)
%
% Computes a loss matrix |S| x n where S is the set of states 
% and n the length of the true state sequence.
%
% true_state_seq -- true sequence of states (of length n) to be learned.
% state_model -- a struct specifying the state transition model (see
%   make_model.m for details)
% PAR -- a struct of parameters specified in setup_hmsvm_training.m and
%   train_hmsvm.m 
% returns the loss matrix used for decoding to obtain the max margin
%   violator
%
% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2007-2008

FP_loss = 1;
FN_loss = 1;


STATES = get_state_set(PAR);

loss = zeros(length(state_model));


%loss([STATES.start STATES.intergenic STATES.stopCodon STATES.stop], ...
%  [STATES.startCodon, STATES.exonic1, STATES.exonic2, STATES.exonic3]) ...
%= FN_loss;

%loss([STATES.startCodon, STATES.exonic1, STATES.exonic2, STATES.exonic3],...
%  [STATES.start STATES.intergenic STATES.stopCodon STATES.stop]) ...
%= FP_loss;


loss(STATES.startCodon, [STATES.intergenic STATES.stopCodon STATES.start STATES.stop ...
  STATES.exonic1, STATES.exonic2, STATES.exonic3]) = FP_loss;

loss(STATES.stopCodon, [STATES.intergenic STATES.startCodon STATES.start STATES.stop ...
  STATES.exonic1, STATES.exonic2, STATES.exonic3]) = FP_loss;


loss(STATES.exonic1, STATES.stopCodon) = FP_loss;
loss(STATES.stopCodon, STATES.exonic1) = FN_loss;

loss(STATES.intergenic, STATES.startCodon) = FN_loss;




for i=1:length(true_state_seq),
  loss_matrix(:,i) = loss(:,true_state_seq(i));
end
