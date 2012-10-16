function mTIM_config
% Global configuration file 

% add 
addpath(sprintf('/mnt/galaxyTools/tools/oqtans/sw_require/utils/'));
addpath(sprintf('/mnt/galaxyTools/tools/oqtans/sw_require/utils/rproc'));
addpath(sprintf('/mnt/galaxyTools/tools/oqtans/mTIM_2/tools/r2007a'));
addpath(sprintf('/mnt/galaxyTools/tools/oqtans/mTIM_2/tools/sotool/src'));
addpath(sprintf('/mnt/galaxyTools/tools/oqtans/mTIM_2/src/model'));
addpath(sprintf('/mnt/galaxyTools/tools/oqtans/mTIM_2/src/utils'));
addpath(sprintf('/mnt/galaxyTools/tools/oqtans/mTIM_2/src/data_preparation'));

engine = 'octave';
if isequal(engine, 'octave'),
  warning('off', 'Octave:precedence-change');
  warning('off', 'Octave:function-name-clash');
  warning('off', '');
  warning('off', 'Octave:num-to-str');
  warning('off', 'Octave:function-name-clash');
  warning('off', 'Octave:divide-by-zero');
  warning('off', 'Octave:future-time-stamp');
  warning('off', 'solve_qp:constraints');
  warning('off', 'Octave:assign-as-truth-value');
else
  warning('off', 'MATLAB:typeaheadBufferOverflow');
end

% make sure no process stops with a debug prompt
global g_ignore_keyboard
g_ignore_keyboard = 1;
