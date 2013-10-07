function init_paths()
% Extends the matlab path to directories needed by the HMSVM toolbox.
%
% see train_hmsvm.m
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008

% locate hmsvm home directory based on the location of this script
this_name = mfilename();
this_path = mfilename('fullpath');
d = this_path(1:end-length(this_name));

% add necessary directories
addpath([d 'losses'])
addpath([d 'native'])
