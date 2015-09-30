function [rDiff_paths] = set_rdiff_paths()
% [difftest_paths] = set_rdiff_paths()

%Initialize rDiff_paths
rDiff_paths='';

%get the path to the source files
BASEPATH_SRC=getenv('RDIFF_SRC_PATH');

%Determine interpreter
if size(ver('Octave'),1)
    INTERPR = 1;
else
    INTERPR = 0;
end


%Add paths to rDiff_paths
if INTERPR
    rDiff_paths=[rDiff_paths ':' [BASEPATH_SRC '/octave/']];
end
rDiff_paths=[rDiff_paths ':' [BASEPATH_SRC '/../mex/']];
rDiff_paths=[rDiff_paths ':' [BASEPATH_SRC '/tests/']];
rDiff_paths=[rDiff_paths ':' [BASEPATH_SRC '/tools/']];
rDiff_paths=[rDiff_paths ':' [BASEPATH_SRC '/tools/read_utils/']];
rDiff_paths=[rDiff_paths ':' [BASEPATH_SRC '/variance/']];
rDiff_paths=[rDiff_paths ':' genpath([BASEPATH_SRC '/locfit/'])];
addpath(rDiff_paths);

