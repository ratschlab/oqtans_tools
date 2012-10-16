function [state_model, A] = make_model_mTIM(PAR, transition_scores)

% [state_model, A] = make_model_mTIM(PAR, [transition_scores])
%
% Combines user specification and automatic completion of the
% state-transition model.
% 
% PAR -- a struct of parameters specified in setup_training.m and
%   train_hmsvm.m
% transition_scores -- optional parameter; if specified, transition scores
%   of A and a_trans will be set accordingly (otherwise initialized to 0)
% returns the fully specified graphical model (state_model) and the 
%   transition matrix (A)
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008

state_model = specify_model(PAR);

if exist('transition_scores', 'var'),
  [state_model, A] = complete_model(state_model, PAR, ...
                                             transition_scores);
else
  [state_model, A] = complete_model(state_model, PAR);
end

% eof