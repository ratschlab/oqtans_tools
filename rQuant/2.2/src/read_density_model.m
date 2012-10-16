function density_model = read_density_model(fname)
% READ_DENSITY_MODEL   Reads density model from text file into matrix.
%
%   density_model = read_density_model(fname)
%
%   -- input --
%   fname:         name of input file 
%
%   -- input --
%   density_model: profiles or intron distances
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2010 Max Planck Society
%


fd = fopen(fname, 'r');

r = 0;
while (1)
  input = fgetl(fd);
  if ~ischar(input), break; end
  items = separate(strtrim(input));
  r = r + 1;
  for c = 1:length(items),
    density_model(r, c) = str2double(items{c});
  end
end
fclose(fd);

if isempty(density_model)
  error('\ncould not read density model\n');
end