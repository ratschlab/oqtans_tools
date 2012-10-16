function main_mll()
% Test the Multiple Loss Learning SVM in its primal form using
% bundle methods.
%
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

% include all paths neccessary to start the demo
addpath('../../src');
addpath('../../src/linesearch');
addpath('../../src/losses');
addpath('../../src/solver/prloqo');

% amount of data points for each class
n1 = 100;
n2 = 100;


% scaling factor for each class
s1 = 10.9;
s2 = 10.9;


% mean values
m1 = [2;1];
m2 = [2;3];

randn('seed',1);

% generate training data
X = s1*randn(2,n1)+repmat(m1,1,n1);
y = ones(1,n1);
X = [X s2*randn(2,n2)+repmat(m2,1,n2)];
y = [y -ones(1,n2)];


% add bias dimension to X 
X = [X; ones(1,length(y))];

% bundle method related parameters
PAR.epsilon = 1e-6;  % (does not change outcome)
PAR.max_num_iter = 250; % bmrm should converge fast
PAR.bmrm.eps = 1e-3; % important: how accurate should the solver be
PAR.bmrm.K = 13; % important: number of cutting planes used
PAR.bmrm.useAggregation = 1; % (do not change)

% parabola line search parameters
PAR.bmrm.lsSteps = 30; % max number of linesearch steps
PAR.bmrm.lsEps = 1e-3; % accuracy
PAR.bmrm.lsTheta = 0.1; % linesearch regularization value


% Loss functions and subgradient method

% Set the function handle for the loss function.
% If possible use native versions of the loss function:
% (see /losses and /native)
% loss_hinge (alt. loss_hinge_native)
% loss_sqrhinge (alt. loss_sqrhinge_native)
% ..
% For multiple loss learning use loss_mll
PAR.loss.fct = @loss_mll; 
% set the corresponding subgradient method 
% For multiple loss learning use loss_mll_sg
PAR.loss.sg_fct = @loss_mll_sg;

% loss specific parameters
PAR.loss.sqrhinge_alpa = 1; % (only sqr hinge) multiplier (sometimes values are too big)

% Params only for multiple loss learning:
% set the loss function handles you want to use
PAR.loss.mll_losses = {@loss_hinge_native, @loss_logistic_native, @loss_sqrhinge_native}
% ..and corresponding subgradients
PAR.loss.mll_losses_sg = {@loss_hinge_sg, @loss_logistic_sg, @loss_sqrhinge_sg}
% ..and the p norm
PAR.loss.mll_p = 2.0;



% train the primal svm with bundle methods
[w, thetas, inds] = train_primal_svm(PAR, 10.01, X, y);
% extract bias term and weight vector
b = w(end);
w = w(1:end-1);

% plot training data and decision plane
figure(1);
subplot(1,2,1);
hold on
plot(X(1,y==1),X(2,y==1),'.b','linewidth',2);
plot(X(1,y==-1),X(2,y==-1),'.r','linewidth',2);

dw = -w(1)/w(2);
db = -b/w(2);
plot([-10 10],[-10*dw+db 10*dw+db],'-c','linewidth',2);

ylim([-5 5]);
xlim([-5 5]);
hold off

subplot(1,2,2);
hold on;
names = {};
for m=1:size(thetas,2),
  name = func2str(PAR.loss.mll_losses{m});
  names{m} = substring(name,5,length(name)-8);
  plot(1:length(y), thetas(inds,m),'color',rand(1,3));
end
title(sprintf('norm(w)=%1.2f',norm(w)));
legend(names);
hold off;
