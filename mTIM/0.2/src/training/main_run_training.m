function main_run_training(CFG)

%
% Runs the whole HM-SVM training including cross-validation
% and model selection.
%
% see HM-SVM toolbox for HM-SVM specifics.
%
% written by Georg Zeller, Pramod Mudrakarta & Gunnar Raetsch, MPI Tuebingen, Germany, 2009-2011



%%% start cross-validation for each parameter configuration
JOB_INFO = [];
for i=1:size(CFG.train_params,1),

  % configure (=set model parameters)
  fprintf('Configuration of model %i...\n\n', i);
  for j=1:length(CFG.param_names),
    if ischar(CFG.train_params{i,j}) ,
      fprintf('  %s = %s\n', CFG.param_names{j}, CFG.train_params{i,j});
    elseif length(CFG.train_params{i,j}) == 1,
      fprintf('  %s = %f\n', CFG.param_names{j}, CFG.train_params{i,j});
    else
      p_str = [];
      for k=1:length(CFG.train_params{i,j}),
        p_str = [p_str sprintf('%f ', CFG.train_params{i,j}(k))];
      end
      fprintf('  %s = %s\n', CFG.param_names{j}, p_str);
    end
    CFG.PAR = setfield(CFG.PAR, CFG.param_names{j}, CFG.train_params{i,j});
  end
  fprintf('\n');
  
  CFG.PAR.model = i;
  if CFG.use_rproc,
    if i==1,
      CFG.rproc_par.identifier = sprintf('%s_m%i', ...
                                         CFG.rproc_par.identifier, i);
    else
      idx = find(CFG.rproc_par.identifier=='_', 1, 'last') - 1;
      CFG.rproc_par.identifier = sprintf('%s_m%i', ...
                                         CFG.rproc_par.identifier(1:idx), i);
    end
  end  
  JOB_INFO = [JOB_INFO, run_cross_validation(CFG, JOB_INFO)];
end

fprintf('Approaching rproc_wait.\n');
if (CFG.use_rproc && CFG.PAR.submit_train),
    % wait for job finished
    fprintf('Waiting for jobs to finish...');
    JOB_INFO = rproc_wait(JOB_INFO, 20, 1, -1);
    fprintf('Done!\n');

    % store the log files in the corresponding training directories
    for i=1:length(JOB_INFO),
        log_fname = JOB_INFO(i).log_fname;
        out_dir = JOB_INFO(i).train_dir;
        fprintf('Copy log file "%s" to target directory "%s"\n',log_fname,out_dir);
        unix(sprintf('cp -f %s %s', log_fname, out_dir)) ; 
    end
end


% eof
