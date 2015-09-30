function mTIM_config
% Global configuration file 

OQTANS = getenv('OQTANS_PATH');
OQTANS_DEP = getenv('OQTANS_DEP_PATH');

% add 
%addpath(sprintf('%s/utils/', OQTANS_DEP));
%addpath(sprintf('%s/utils/rproc', OQTANS_DEP));
addpath(sprintf('%s/mTIM/0.2/tools', OQTANS));
addpath(sprintf('%s/mTIM/0.2/tools/r2007a', OQTANS));
addpath(sprintf('%s/mTIM/0.2/tools/sotool', OQTANS));
addpath(sprintf('%s/mTIM/0.2/tools/sotool/losses', OQTANS));
addpath(sprintf('%s/mTIM/0.2/tools/sotool/native', OQTANS));
addpath(sprintf('%s/mTIM/0.2/src/model', OQTANS));
addpath(sprintf('%s/mTIM/0.2/src/evaluation', OQTANS));
addpath(sprintf('%s/mTIM/0.2/src/training', OQTANS));
addpath(sprintf('%s/mTIM/0.2/src/prediction', OQTANS));
addpath(sprintf('%s/mTIM/0.2/src/utils', OQTANS));
addpath(sprintf('%s/mTIM/0.2/src/data_preparation', OQTANS))

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
