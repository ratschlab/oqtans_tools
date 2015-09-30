function JOB_INFO = run_cross_validation(CFG, JOB_INFO)
%
% run_cross_validation(CFG, JOB_INFO)
%
% Starts cross-validation.
%
% see HM-SVM toolbox for HM-SVM specifics.
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2011


fprintf('  Starting %i-fold cross_validation...\n\n', CFG.num_xval_folds);
disp(CFG);

for f=1:CFG.num_xval_folds,
  % input data file
  CFG.PAR.data_file = [CFG.xval_dirs{f} '/train_data']; 
  % output directory 
  CFG.PAR.out_dir = CFG.xval_train_dirs{CFG.PAR.model,f};
  
  fprintf('    Starting training for fold %i...\n\n', f);
  disp(CFG.PAR);
  
  % partition exon intensities into expression bins 
  % (i.e. assign discrete expression levels)
  clear signal
  if (~exist(CFG.PAR.data_file)),
    load([CFG.PAR.data_file '.mat'], 'signal', 'label', 'exm_id_intervals');
  else
    load(CFG.PAR.data_file, 'signal', 'label', 'exm_id_intervals');
  end

  if ~exist('signal', 'var'),
    fn_data = CFG.PAR.data_file;
    if isequal(fn_data(end-3:end), '.mat'),
      fn_data = fn_data(1:end-4);
    end
    fn_data = [fn_data '_signal'];
    if(~exist(fn_data)),
        signal = load_struct([fn_data '.mat'], 'signal');
    else
        signal = load_struct(fn_data, 'signal');
    end
  end
  % compute expression bins for th expression-specific submodels based on
  % exon read coverage (averaged over + and - strand)
  FEATS = get_feature_set_mTIM();
  EC = strmatch('exon_cover', FEATS, 'exact');
  bins = discretize_expression(signal(EC,:), label, exm_id_intervals, CFG.PAR.num_levels);
  assert(size(bins,1) == CFG.PAR.num_levels);
  CFG.PAR.expression_bins = bins;
  
  % filter potential training examples
  n = CFG.PAR.num_train_exm;
  if (n==-1),
%      warning('Using only minimum training set filter!');
%      CFG.PAR.train_exms = min_filter_training_data(label, signal, exm_id_intervals, CFG);
      CFG.PAR.train_exms = filter_training_data(label, signal, exm_id_intervals, CFG);
  else
      CFG.PAR.train_exms = filter_training_data(label, signal, exm_id_intervals, CFG);
  end


  % set the number of features to be used for leraning
  CFG.PAR.num_features = size(signal, 1);
  clear signal 


  % subselect sequences for training
  if (n==-1),
    % take all examples for training
    fprintf('Choosing from %i filtered training examples\n', length(CFG.PAR.train_exms));
    fprintf('Using ALL training examples\n');
    fprintf('Disabling validation examples\n', CFG.num_vald_exm);
    CFG.PAR.train_exms = CFG.PAR.train_exms;
    CFG.PAR.vald_exms = [];
    CFG.num_vald_exm = 0;
  else
    % take some trainig sequences, remove them from training set 
    % and instead use them for ad hoc accuracy estimation during training
    holdout = randperm(length(CFG.PAR.train_exms));
    holdout = holdout(1:CFG.num_vald_exm);
    CFG.PAR.vald_exms = CFG.PAR.train_exms(holdout);
    CFG.PAR.train_exms(holdout) = [];
    fprintf('Choosing from %i filtered training examples\n', length(CFG.PAR.train_exms));
    fprintf('Using %i training examples\n', n);
    fprintf('Using %i validation examples\n', CFG.num_vald_exm);
    r = randperm(length(CFG.PAR.train_exms));
    CFG.PAR.train_exms = CFG.PAR.train_exms(r(1:n));
  end
  assert(isempty(intersect(CFG.PAR.train_exms, CFG.PAR.vald_exms)));

  exm_ids = sort([CFG.PAR.train_exms; CFG.PAR.vald_exms]);
  assert(all(ismember(exm_ids, [CFG.PAR.train_exms; CFG.PAR.vald_exms])));
  % re-save training data which will actually be used
  CFG.PAR.data_file = save_reduced_data(CFG.PAR.data_file, exm_ids, CFG);
 
  bak = CFG.PAR.num_train_exm;
  CFG.PAR.num_train_exm = length(CFG.PAR.train_exms);

  %if CFG.grid_use,
    % don't let the matgrid toolbox collect the results
  %  JOB.id = mgsub({}, 'train_hmsvm', {CFG.PAR}, 'qsub_opts', sprintf('-l h_vmem=%iG', CFG.grid_memreq));
  %  JOB.train_dir = CFG.PAR.out_dir;
  %  JOB_INFO = [JOB_INFO, JOB];    % add information about the training output directory
  %  fprintf('    Submitted job %i\n\n', length(JOB_INFO));
  %else
    % train sequentially for different parameter combinations
    train_hmsvm(CFG.PAR);
  %end

  CFG.PAR.num_train_exm = bak;
end
% eof
