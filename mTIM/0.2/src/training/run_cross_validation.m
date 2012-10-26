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
    load(CFG.PAR.data_file, '-mat', 'signal', 'label', 'exm_id_intervals');
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
  CFG.PAR.train_exms = filter_training_data(label, signal, exm_id_intervals, CFG);

  % set the number of features to be used for leraning
  CFG.PAR.num_features = size(signal, 1);
  clear signal 

  % take some trainig sequences, remove them from training set 
  % and instead use them for ad hoc accuracy estimation during training
  holdout = randperm(length(CFG.PAR.train_exms));
  holdout = holdout(1:CFG.num_vald_exm);
  CFG.PAR.vald_exms = CFG.PAR.train_exms(holdout);
  CFG.PAR.train_exms(holdout) = [];
  assert(isempty(intersect(CFG.PAR.train_exms, CFG.PAR.vald_exms)));

  % subselect sequences for training
  n = CFG.PAR.num_train_exm;
  fprintf('Choosing from %i filtered training examples\n', length(CFG.PAR.train_exms));
  fprintf('Using %i training examples\n', n);
  fprintf('Using %i validation examples\n', CFG.num_vald_exm);
  r = randperm(length(CFG.PAR.train_exms));
  CFG.PAR.train_exms = CFG.PAR.train_exms(r(1:n));
  exm_ids = sort([CFG.PAR.train_exms; CFG.PAR.vald_exms]);
  assert(all(ismember(exm_ids, [CFG.PAR.train_exms; CFG.PAR.vald_exms])));
  % re-save training data which will actually be used
  CFG.PAR.data_file = save_reduced_data(CFG.PAR.data_file, exm_ids, CFG);
  
  if CFG.use_rproc && CFG.PAR.submit_train,
    % submit cluster job
    if f==1,
      CFG.rproc_par.identifier = sprintf('%s_f%i_', ...
                                         CFG.rproc_par.identifier, f);
    else
      idx = find(CFG.rproc_par.identifier=='_', 2, 'last') - 1;
      CFG.rproc_par.identifier = sprintf('%s_f%i_', ...
                                         CFG.rproc_par.identifier(1:idx(1)), f);
    end

    
    JOB =  rproc('train_hmsvm', CFG.PAR, CFG.rproc_memreq, CFG.rproc_par, CFG.rproc_time);
    JOB.train_dir = CFG.PAR.out_dir;
    JOB_INFO = [JOB_INFO, JOB];    % add information about the training output directory
    fprintf('    Submitted job %i\n\n', length(JOB_INFO));
  else
    
    % train sequentially for different parameter combinations
    train_hmsvm(CFG.PAR);
  end
end

% eof
