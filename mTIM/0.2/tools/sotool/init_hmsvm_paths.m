function init_hmsvm_paths
% init only the paths corresponding to the hmsvm
% NOT the bmrm training paths
% (needed e.g. parallel prediction when calling
% decode_viterbi directly).

% add paths
this_name = mfilename();
this_path = mfilename('fullpath');
idx = find(this_path=='/', 1, 'last');
src_dir = this_path(1:idx);
fprintf('Hmsvm-source directory is %s.\n',src_dir);

% add sotool-app paths
addpath([src_dir '/native']);


