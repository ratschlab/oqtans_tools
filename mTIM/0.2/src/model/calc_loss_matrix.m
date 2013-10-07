function loss_matrix = calc_loss_matrix(true_state_seq, state_model, PAR)

% loss_matrix = calc_loss_matrix(true_state_seq, state_model)
%
% Computes a loss matrix |S| x n where S is the set of states 
% and n the length of the true state sequence.
%
% true_state_seq -- true sequence of states (of length n) to be learned.
% state_model -- a struct specifying the state transition model (see
%   make_model_mTIM.m for details)
% PAR -- a struct of parameters specified in setup_training.m and
%   train_hmsvm.m 
% returns the loss matrix used for decoding to obtain the max margin
%   violator, row index corresponds to predicted state, column index to
%   true state
%
% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2008-2011

fp_loss = 0.6;
fn_loss = 0.4; % TODO: was 0.4 for the good version pred2011-04-28_14h56

% 2nd try
% TODO: 10.10.2012
% 1. set the false positive penelization lower to increase sensitivity
% 2. wrong_level_loss is 0.05 for all levels (=increase specificity for low levels)
% 3. increase intron to ige penelization in order to decrease amount of splitted genes
% 4. increase exon to ige penelization in order to decrease amount of splitted genes
fp_loss = 0.5;
fn_loss = 0.5;


if isfield(PAR, 'exon_intron_loss'),    
  exon_intron_loss = PAR.exon_intron_loss;
else
  exon_intron_loss      = 0.1;
end
if isfield(PAR, 'strand_loss'),    
  strand_loss = PAR.strand_loss;
else
  strand_loss           = 0.1;
end
if isfield(PAR, 'level_loss') && PAR.level_loss == 1,    
  level_loss = log([2:PAR.num_levels+1]);
  level_loss = level_loss / mean(level_loss);
else
  level_loss = ones(1,PAR.num_levels);
end
% TODO: 10.10.2012
wrong_level_loss        = 0.01;
wrong_level_loss        = 0.01;
intron_boundary_loss    = 10; % TODO: was 5 for the good version pred2011-04-28_14h56

STATES = get_state_set_mTIM(PAR);

fn = fieldnames(STATES);
efw_idx = strmatch('EFW', fn);
eiw_idx = strmatch('EIW', fn);
elw_idx = strmatch('ELW', fn);
exw_idx = sort([strmatch('EFW', fn); strmatch('EIW', fn); strmatch('ELW', fn)]);
assert(length(exw_idx) == 3*PAR.num_levels);
efc_idx = strmatch('EFC', fn);
eic_idx = strmatch('EIC', fn);
elc_idx = strmatch('ELC', fn);
exc_idx = sort([strmatch('EFC', fn); strmatch('EIC', fn); strmatch('ELC', fn)]);
assert(length(exc_idx) == 3*PAR.num_levels);

ifw_idx = strmatch('IFW', fn);
iiw_idx = strmatch('IIW', fn);
ilw_idx = strmatch('ILW', fn);
inw_idx = sort([strmatch('IFW', fn); strmatch('IIW', fn); strmatch('ILW', fn)]);
assert(length(inw_idx) == 3*PAR.num_levels);

ifc_idx = strmatch('IFC', fn);
iic_idx = strmatch('IIC', fn);
ilc_idx = strmatch('ILC', fn);
inc_idx = sort([strmatch('IFC', fn); strmatch('IIC', fn); strmatch('ILC', fn)]);
assert(length(inc_idx) == 3*PAR.num_levels);

loss = zeros(length(state_model));

% Have a loss for exon and intron predictions on the wrong strand
loss([exw_idx' inw_idx'], exc_idx) = strand_loss;
loss([exw_idx' inw_idx'], inc_idx) = strand_loss;
loss([exc_idx' inc_idx'], exw_idx) = strand_loss;
loss([exc_idx' inc_idx'], inw_idx) = strand_loss; 

% Have a loss for false positive exon and intron predictions
% TODO: 10.10.2012
loss([exw_idx' exc_idx'], STATES.IGE) = fp_loss;


% TODO: 10.10.2012
%loss([inw_idx' inc_idx'], STATES.IGE) = 0.8*fp_loss; % TODO: was 1 for the good version pred2011-04-28_14h56
loss([inw_idx' inc_idx'], STATES.IGE) = 1.0*fp_loss; % TODO: was 1 for the good version pred2011-04-28_14h56
 
% Have the loss for false negative predictions (exons & introns) weighted by
% expression level, i.e. penalize low expressed exons that go undetected
% less than highly expressed ones
for i=1:PAR.num_levels,
  loss(STATES.IGE, eiw_idx(i)) = level_loss(i)*fn_loss;
  loss(STATES.IGE, eic_idx(i)) = level_loss(i)*fn_loss;
  loss(STATES.IGE, iiw_idx(i)) = level_loss(i)*fn_loss;
  loss(STATES.IGE, iic_idx(i)) = level_loss(i)*fn_loss;
end

% Have a loss for exon intron / confusions
loss(exw_idx', inw_idx') = exon_intron_loss;
loss(exc_idx', inc_idx') = exon_intron_loss;
loss(inw_idx', exw_idx') = exon_intron_loss;
loss(inc_idx', exc_idx') = exon_intron_loss;

% Have a loss for false-positive intron boundary states
loss([ifw_idx'], [STATES.IGE iiw_idx' ilw_idx' exw_idx']) = intron_boundary_loss;
loss([ifw_idx'],            [iic_idx' ilc_idx' exc_idx']) = intron_boundary_loss;
loss([ifc_idx'], [STATES.IGE iic_idx' ilc_idx' exc_idx']) = intron_boundary_loss;
loss([ifc_idx'],            [iiw_idx' ilw_idx' exw_idx']) = intron_boundary_loss;
loss([ilw_idx'], [STATES.IGE iiw_idx' ifw_idx' exw_idx']) = intron_boundary_loss;
loss([ilw_idx'],            [iic_idx' ifc_idx' exc_idx']) = intron_boundary_loss;
loss([ilc_idx'], [STATES.IGE iic_idx' ifc_idx' exc_idx']) = intron_boundary_loss;
loss([ilc_idx'],            [iiw_idx' ifw_idx' exw_idx']) = intron_boundary_loss;

 % TODO: was 0 (non-existent) for the good version pred2011-04-28_14h56
% Have a loss for false-negative intron boundary states

% TODO: 10.10.2012
FN_INTRON_LOSS = 0.2;
FN_INTRON_LOSS = 0.2;
loss([STATES.IGE iiw_idx' ilw_idx' exw_idx'], [ifw_idx']) = FN_INTRON_LOSS*intron_boundary_loss;
loss([iic_idx' ilc_idx' exc_idx'],            [ifw_idx']) = FN_INTRON_LOSS*intron_boundary_loss;
loss([STATES.IGE iic_idx' ilc_idx' exc_idx'], [ifc_idx']) = FN_INTRON_LOSS*intron_boundary_loss;
loss([iiw_idx' ilw_idx' exw_idx'],            [ifc_idx']) = FN_INTRON_LOSS*intron_boundary_loss;
loss([STATES.IGE iiw_idx' ifw_idx' exw_idx'], [ilw_idx']) = FN_INTRON_LOSS*intron_boundary_loss;
loss([iic_idx' ifc_idx' exc_idx'],            [ilw_idx']) = FN_INTRON_LOSS*intron_boundary_loss;
loss([STATES.IGE iic_idx' ifc_idx' exc_idx'], [ilc_idx']) = FN_INTRON_LOSS*intron_boundary_loss;
loss([iiw_idx' ifw_idx' exw_idx'],            [ilc_idx']) = FN_INTRON_LOSS*intron_boundary_loss;


% Also penalize if expression levels are confused for exon states
% in a manner increasing with level difference
level_lm = zeros(PAR.num_levels);

% TODO: 10.10.2012
%for i=1:PAR.num_levels-1,
i = 1;
for n=1:PAR.num_levels-1,
  i = n;
  level_lm = level_lm + diag(repmat(i*wrong_level_loss,1,PAR.num_levels-n), +n);
  level_lm = level_lm + diag(repmat(i*wrong_level_loss,1,PAR.num_levels-n), -n);
end
loss(eiw_idx, eiw_idx) = level_lm;
loss(eic_idx, eic_idx) = level_lm;


assert(all(diag(loss))==0);

%imagesc(loss);
%xlabel('true state');
%ylabel('predicted state');
%keyboard

loss_matrix = compute_loss_matrix(loss, true_state_seq);

% simple hamming loss:
% of course, performance drops, but not as much as expected...
%
%loss_matrix = ones(size(loss,1),length(true_state_seq));
%for i=1:length(true_state_seq),
%  loss_matrix(true_state_seq(i),i) = 0.0;
%end


if PAR.extra_checks,
  l = 0;
  for i=1:length(true_state_seq),
    l = l + loss_matrix(true_state_seq(i),i);
  end
  assert(l==0);
end

% eof
