function CFG = create_dirs(CFG,create_prep_dirs,create_exp_dirs)
% Creates all directories for a single experiment.
%
%

fprintf('\nCreating directories for experiment %s...\n\n',CFG.exp_name);

if (~exist(CFG.out_train_dir,'dir') && create_exp_dirs), 
    CFG.out_train_dir
    mkdir(CFG.out_train_dir); 
end;
if (~exist(CFG.out_pred_dir,'dir') && create_exp_dirs), 
    mkdir(CFG.out_pred_dir); 
end;

% create model directories for training
CFG.model_dirs = {};
for i=1:size(CFG.train_params,1),
  d = sprintf('%s/model%i/', CFG.out_train_dir, i);
  if (~exist(d,'dir') && create_exp_dirs), 
      fprintf('Creating model-directory %s.\n',d);
      mkdir(d); 
  end
  CFG.model_dirs{i} = d;  
end



% Create the cross-validation directories (for data_preparation and
% for each model of the training)
CFG.xval_dirs = {};
CFG.xval_train_dirs = {};
for fold=1:CFG.num_xval_folds,
  d = sprintf('%s/xval_fold%i/', CFG.out_dir, fold);
  if (~exist(d,'dir') && create_prep_dirs), 
      fprintf('Creating xval-directory %s.\n',d);
      mkdir(d); 
  end
  CFG.xval_dirs{fold} = d;  

  % for each model
  for model=1:length(CFG.model_dirs),
    d = sprintf('%s/xval_fold%i/', CFG.model_dirs{model}, fold);
    if (~exist(d,'dir') && create_exp_dirs), 
      fprintf('Creating training xval-directory %s.\n',d);
      mkdir(d); 
    end
    CFG.xval_train_dirs{model,fold} = d;  
  end
end


