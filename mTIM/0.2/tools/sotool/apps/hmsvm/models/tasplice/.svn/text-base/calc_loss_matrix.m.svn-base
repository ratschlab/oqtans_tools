function loss_matrix = calc_loss_matrix(true_state_seq, state_model, PAR)
% loss_matrix = calc_loss_matrix(true_state_seq, state_model, PAR)
% compute loss matrix |S| x n where S is the set of states 
% and n the length of the true state sequence

% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany

FP_segment_loss = 1;
FN_segment_loss = 1;

STATES = get_state_set(PAR);

loss = zeros(length(state_model));
loss([STATES.ino STATES.ige], STATES.exo) = FN_segment_loss;
loss(STATES.exo, [STATES.ino STATES.ige]) = FP_segment_loss;

%loss([STATES.ino STATES.ige STATES.exo], ...
%     [STATES.ino_ss STATES.ige_ss STATES.exo_ss STATES.trans_start STATES.trans_end STATES.don STATES.acc]) ...
%    = ;
%loss(STATES.exo, [STATES.ino STATES.ige]) ...
%    = FP_loss;

for i=1:length(true_state_seq),
  loss_matrix(:,i) = loss(:,true_state_seq(i));
end
