function main_run_prediction(CFG)
% Main function for making whole genome predictions with mTIM.
%
% see HM-SVM toolbox for HM-SVM specifics.
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009


num_models = size(CFG.xval_train_dirs,1);
fprintf('There are %i models present.\n',num_models);

% for any model
for m=1:num_models,

    % discard previous predictions
    genes = [];
    vald_genes = [];

    % for predicting the wohle genome we only need one predictor (and thus
    % pick that from the first cross-validation fold) 
    num_predictors_needed = CFG.num_xval_folds;
    if (CFG.predict_genome),
      num_predictors_needed = 1;
      warning(['Only using predictor 1 to generate whole-genome predictions\nIf ' ...
               'TEST predictions are to be made, set CFG.predict_genome to 0!']);
    end
    for f=1:num_predictors_needed,
      predictor_dir = CFG.xval_train_dirs{m,f};
      
      % predict test data
      CFG.predict_vald = 0;
      genes = [genes predict_fold(predictor_dir, CFG)];

      if (CFG.use_vald_data),
          % predict validation data
          CFG.predict_vald = 1;
          vald_genes = [vald_genes predict_fold(predictor_dir, CFG)];
      end
    end

    if ~isempty(genes),
        % sort genes (necessary because test chunks are not sorted)
        tmp = [[genes.chr_num]', [genes.start]', [genes.stop]'];
        [tmp, idx] = sortrows(tmp);
        genes = genes(idx);
    end

    % save genes
    fn = [CFG.model_dirs{m} CFG.pred_filename];
    fprintf('Saving test chunk predictions to %s.\n',fn);
    save(fn, 'genes', 'CFG');
    
    % evaluate prediction
    CFG.eval_vald=0;
    CFG.opt_model=m;
    main_run_evaluation(CFG);


    % save validation genes
    if (CFG.use_vald_data),
        genes = vald_genes;
        if ~isempty(genes),
            % sort genes (necessary because vald chunks are not sorted)
            tmp = [[genes.chr_num]', [genes.start]', [genes.stop]'];
            [tmp, idx] = sortrows(tmp);
            genes = genes(idx);
        end

        fn = [CFG.model_dirs{m} CFG.vald_filename];
        fprintf('Saving vald chunk predictions to %s.\n',fn);
        save(fn, 'genes', 'CFG');
    
        % evaluate validation
        CFG.eval_vald=1;
        CFG.opt_model=m;
        main_run_evaluation(CFG);
   end

end

% TODO: for a proper model selection it is
% now the time to select the best prediction
% based on validation for each measure (intron-f1,gene.prec,...)
warning('TODO: model selection: select best predictor for each measure.');

% eof



