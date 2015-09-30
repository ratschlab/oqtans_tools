function main_run_training(CFG)
%
% Runs the whole HM-SVM training including cross-validation
% and model selection.
%
% see HM-SVM toolbox for HM-SVM specifics.
%
% written by Georg Zeller, Pramod Mudrakarta & Gunnar Raetsch, MPI Tuebingen, Germany, 2009-2011
% adapted by Nico Goernitz, Berlin Institute of Technology, 2012
%    - migration from rproc to matgrid toolbox


%%% start cross-validation for each parameter configuration
JOB_INFO = [];
for i=1:size(CFG.train_params,1),

  % Configuration of the parameter struct (=set model parameters)
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

  % run the 
  CFG.PAR.model = i;
  JOB_INFO = [JOB_INFO, run_cross_validation(CFG, [])];
end


fprintf('Approaching parallelization idle loop.\n');
if (CFG.grid_use),
    % wait for job finished
    fprintf('Waiting for jobs to finish...');
    while true,
        if mgshow([JOB_INFO.id],1,30), break; end;
    end
    fprintf('Done!\n');

    % store the log files in the corresponding training directories
    for i=1:length(JOB_INFO),
        log_dir = mgjobdir(JOB_INFO(i).id);
        jobid = JOB_INFO(i).id;
        gridid = mgsungridid(jobid);
        log_fname = sprintf('%s/mgjob-%i.sh.o%i', log_dir, jobid, gridid);
        out_dir = JOB_INFO(i).train_dir;
        
        fprintf('Copy log file "%s" to target directory "%s"\n',log_fname,out_dir);
        unix(sprintf('cp -f %s %s', log_fname, out_dir)) ; 
    end
end


% eof
