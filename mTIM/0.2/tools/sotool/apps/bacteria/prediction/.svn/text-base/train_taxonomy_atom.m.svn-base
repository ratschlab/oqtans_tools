function train_taxonomy_atom(PAR)
%
%
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

addpath('../model');

addpath('../../hmsvm');
addpath('../../hmsvm/native');

addpath('../../../src');
addpath('../../../src/linesearch');
addpath('../../../src/losses');
addpath('../../../src/solver/prloqo');

model = PAR.model;
trainset = PAR.trainset;
PAR.model = [];
PAR.trainset = [];
[information] = train_primal_sosvm(PAR, model, trainset);
% eof

