function evaluate_result(dirname)
% Evaluate the final.mat result in 'dirname'. Stores another
% file 'evaluation.mat' in the same directory which can be
% use as multi-task file.
%
% written by Nico Goernitz, MPI Tuebingen, TU Berlin, Germany, 2011

% load neccessary directories
addpath('model');
addpath('../hmsvm');
addpath('../hmsvm/native');

%dirname='../../out/CTEST';
fprintf('Evaluate result in %s.\n',dirname);

% load final result file
assert(exist(sprintf('%s/final.mat',dirname),'file')>=1);
load([dirname '/final.mat']);

% number of repetitions and parameters
% can be extracted from final.mat result variable
assert(exist('result','var')>=1);
num_orgs = length(PAR.organism);
num_reps = size(result,1);
num_params = size(result,2);
fprintf('%i organisms with %i repetitions and %i parameters.\n',num_orgs,num_reps,num_params);

% check amounts of test data for each organism
num_org_test = [];
org_name = {};
org_test = {};
org_vald = {};
for i=1:num_orgs,
  org_name{i} = PAR.organism{i}.name;
  org_test{i} = PAR.organism{i}.test_exms;
  org_vald{i} = PAR.organism{i}.vald_exms;
  num_org_test(i) = length(org_test{i});
  fprintf('Organism %s has %i test exm.\n', PAR.organism{i}.name, num_org_test(i));
end

% idea 1: reweight the test results for each organism because
% of different test set sizes
% weight for each organism
weight = num_org_test ./ max(num_org_test);
weight = 1./weight;  

% idea 2: only use min(num_org_test) examples for each
% organism
num_test = min(num_org_test);
fprintf('Use %i test examples for each organism.\n',num_test);

% choose the best classifier on base of the complete validation data set
% (number of vald data is equal for every organism).
out = zeros(num_params,4);
inds = zeros(num_params,1);
for n=1:num_reps,
  load([dirname sprintf('/data_%i.mat',n)]);
  for i=1:num_params,
    PAR = result{n,i}.PAR;
    model = result{n,i}.model;
    % same amount of validation examples for each organism
    res = eval_set(PAR, PAR.vald_exms, model, exm_id_intervals, signal, label);
    out(i,:) = out(i,:) + [res.sensitivity res.specificity res.precision res.f_score];
    
    % sign all parameters which have mtl_b=1.0 (==no mtl)
    % with mtl enabled and only one organism (therefore this
    % is a leaf node)
    inds(i) = num_orgs==1 && PAR.mtl.mtl_enable && abs(PAR.mtl.mtl_b-1.0)<=1e-8;
  end
end

% get best model
out = out/num_reps;
fprintf('[sensitivity specificity precision f_score mtl_b==1.0]\n');
disp([out, inds]);

% best score for all parameters
[best_score idx] = max(out(:,4));
fprintf('Best f_score(%1.4f) with C=%3.2f mtl_b(enabled)=%1.4f(%i)\n', ...
  best_score, result{1,idx}.PAR.C_small, result{1,idx}.PAR.mtl.mtl_b, result{1,idx}.PAR.mtl.mtl_enable);

% what is the result for the single organism without MTL?
res_single = [];
inds = find(logical(inds));
[best_score_single idx_single] = max(out(inds,4));
idx_single = inds(idx_single);

if (~isempty(idx_single)),
  fprintf('Best SINGLE f_score(%1.4f) with C=%3.2f mtl_b(enabled)=%1.4f(%i)\n', ...
    best_score_single, result{1,idx_single}.PAR.C_small, result{1,idx_single}.PAR.mtl.mtl_b, PAR.mtl.mtl_enable);

  res_single.name = org_name{1};
  res_single.result = eval_set(result{end,idx_single}.PAR, org_test{1}, ...
    result{end,idx_single}.model, exm_id_intervals, signal, label);
  % some feedback
  fprintf('SINGLE Organism %s has %i test examples and achieved performance(f_score): %1.4f\n',...
    org_name{1}, length(org_test{1}), res_single.result.f_score);

end

w = result{end,idx}.w;
PAR = result{end,idx}.PAR;
model = result{end,idx}.model;

% for each organism load 'num_test' test ids
ids = [];
res = {};
for i=1:num_orgs,
  foo = org_test{i};
  res{i}.name = org_name{i};
  res{i}.result = eval_set(PAR, foo, model, exm_id_intervals, signal, label);
  % some feedback
  fprintf('Organism %s has %i test examples and achieved performance(f_score): %1.4f\n',...
    org_name{i}, length(foo), res{i}.result.f_score);
  
  % collect ids for overall testing
  ids = [ids foo(randperm(num_test))];
end

% overall test (only if more than 1 organism)
res_all = [];
if (num_orgs>1),
  res_all.name = 'overall performance';
  res_all.result = eval_set(PAR, ids, model, exm_id_intervals, signal, label);
  fprintf('Overall performance (with %i test examples each organism, f_score): %1.4f\n',num_test, res_all.result.f_score);
end

% save the results
fprintf('Saving result to evaluation.mat...');
save([dirname '/evaluation.mat'],'w','PAR','model','res','res_all','res_single',...
  'best_score','idx','out','splitting');
fprintf('DONE!\n\n');






function res = eval_set(PAR, ids, model, exm_id_intervals, signal, label)
% Evaluate the example with id 'ids'.

% get valid labels
LABELS = get_label_set();

tp=0; tn=0; fp=0; fn=0;
for j=1:length(ids),
  % get the index of the current example id
  idx = find(exm_id_intervals(:,1) == ids(j));
  % start and stop position in the sequence
  % and observation sequence and corresponding label sequence
  start=exm_id_intervals(idx,2);
  stop=exm_id_intervals(idx,3);
  obs_seq = signal(start:stop);
  label_seq = label(start:stop);    

  % prediction
  [pred_path] = decode_viterbi(obs_seq, model.transition_scores, model.score_plifs, PAR);
  
  % evaluate the nukleotide exonic performance
  [f1,f2,f3,f4] = eval_performance(label_seq,pred_path.label_seq,LABELS.exonic);
  tp=tp+f1; tn=tn+f2; fp=fp+f3; fn=fn+f4;
end
% get sensitivity, specificity, f_score, precision
res = eval_confusion_matrix(tp,tn,fp,fn);
