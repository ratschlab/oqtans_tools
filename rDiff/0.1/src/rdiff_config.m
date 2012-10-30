function rdiff_config
% RDIFF_CONFIG   Sets a few global variables with system dependent paths.

global RDIFF_PATH RDIFF_SRC_PATH
global RDIFF_VERSION

% interpreter paths
global INTERPRETER MATLAB_BIN_PATH OCTAVE_BIN_PATH
% SAMTools path
global SAMTOOLS_DIR

% configuration (adapt to the user's configuration)
RDIFF_VERSION = getenv('RDIFF_VERSION');
RDIFF_PATH = getenv('RDIFF_PATH');
RDIFF_SRC_PATH = getenv('RDIFF_SRC_PATH');
INTERPRETER = getenv('INTERPRETER');
MATLAB_BIN_PATH = getenv('MATLAB_BIN_PATH');
OCTAVE_BIN_PATH = getenv('OCTAVE_BIN_PATH');
SAMTOOLS_DIR = getenv('SAMTOOLS_DIR');
OQTANS_DEP_PATH = getenv('OQTANS_DEP_PATH');

% switch off a few expected warnings
addpath(sprintf('%s/utils', OQTANS_DEP_PATH));
addpath(sprintf('%s/utils/rproc', OQTANS_DEP_PATH));
addpath(sprintf('%s/mex', RDIFF_PATH));

engine = determine_engine();
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
