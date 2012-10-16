function rquant_gendata(exm_name, fasta_fname, gff3_fname, alignment_fname, info_fname, profiles_fname)
% RQUANT_GENDATA Loads rQuant example data.
%
%   rquant_gendata(exm_name, fasta_fname, gff3_fname, alignment_fname, info_fname, profiles_fname)
%
%   -- input --
%   exm_name:           name of example data
%   fasta_fname:        name of output fasta file
%   gff3_fname:         name of output gff3 file
%   alignment_fname:    name of output SAM file
%   info_fname:         name of output information file
%   profiles_fname:     name of output file that stores profiles
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2011 Max Planck Society
%


global RQUANT_PATH 

unix(sprintf('cp -p %s/examples/data/%s.fasta %s', RQUANT_PATH, exm_name, fasta_fname));
unix(sprintf('cp -p %s/examples/data/%s.gff3 %s', RQUANT_PATH, exm_name, gff3_fname));
unix(sprintf('cp -p %s/examples/data/%s.sam %s', RQUANT_PATH, exm_name, alignment_fname));
unix(sprintf('cp -p %s/examples/data/profiles.txt %s', RQUANT_PATH, profiles_fname));
unix(sprintf('cp -p %s/examples/data/%s.info %s', RQUANT_PATH, exm_name, info_fname));

return


function ret = keyboard_allowed()
% ret = keyboard_allowed()
%
% returns 1, if a keyboard command would be allowed

global g_ignore_keyboard
global THIS_IS_A_RPROC_PROCESS  

if isequal(g_ignore_keyboard, 1),
  ret = 0;
  return;
end

if isequal(THIS_IS_A_RPROC_PROCESS, 1),
  ret = 0;
  return;
end

ret = 1;

return;


