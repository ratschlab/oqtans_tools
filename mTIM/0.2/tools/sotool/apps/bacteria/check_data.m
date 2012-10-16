function check_data(dirs)
%
%
% adapted by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011



for i=1:length(dirs),
  % partition data for cross-validation
  data_file = [dirs{i} '/data.mat'];

  load(data_file);
  all_exm_ids = unique(exm_id);
  
  fprintf('%s has %i examples.\n',dirs{i}, length(all_exm_ids)); 
end
% eof

