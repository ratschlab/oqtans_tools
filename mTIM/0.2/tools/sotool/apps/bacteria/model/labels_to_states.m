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
% adapted by Nico Goernitz, MPI Tuebingen, TU Berlin, Germany, 2011
problems = [];

LABELS = get_label_set();
STATES = get_state_set(PAR);

% start with complete intergenic labeled sequence
state_seq = repmat(STATES.intergenic, 1, length(label_seq));

% add start/stop codons at specific loci
state_seq(label_seq==LABELS.startCodon) = STATES.startCodon;
state_seq(label_seq==LABELS.stopCodon) = STATES.stopCodon;

% find all start and end codons 
sinds = find(label_seq==LABELS.startCodon);
tinds = find(label_seq==LABELS.stopCodon);
% for each start codon there should be one stop codon
assert(length(sinds)==length(tinds));
% there should be only one gene

%assert(length(sinds)==1);

% we want to include more than one gene for
% each example
sinds = sort(sinds);
tinds = sort(tinds);
for i=1:length(sinds),
  % also handle 'empty' sequences (training examples
  % which have no gene at all
  idx = sinds(i)+1;
  eidx = tinds(i);
  
  % assume that the start codon comes first
  assert(isempty(idx) | idx<eidx);
  % mod-counter for exonic states
  cnt = 1; % after start codon there should be the exonic2 state
  for j=idx:eidx-1;
    state = STATES.exonic1;
    if (cnt==1), state=STATES.exonic2; end;
    if (cnt==2), state=STATES.exonic3; end;    
    state_seq(j) = state;
    cnt = mod(cnt+1,3);
  end
  % check if the state before stop-codon is the exonic3 state and
  % the next is intergenic
  assert(isempty(idx) | state_seq(eidx)==STATES.stopCodon);
  assert(isempty(idx) | state_seq(eidx-1)==STATES.exonic3);
  assert(isempty(idx) | state_seq(eidx+1)==STATES.intergenic);
  % check for the sequence: intergenic startCodon exonic2
  assert(isempty(idx) | (state_seq(idx-1)==STATES.startCodon ...
    & state_seq(idx)==STATES.exonic2 ...
    & state_seq(idx-2)==STATES.intergenic));
end

% set start and stop sequence states
state_seq(1) = STATES.start;
state_seq(end) = STATES.stop;
