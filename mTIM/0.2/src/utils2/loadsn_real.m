function [mat] = loadsn(FILENAME)
% LOADSN load a matrix in SN file format
%	LOADSN (FILENAME) loads an SN matrix (binary or ascii) from
%	the file with name FILENAME

fid = fopen(FILENAME);

fprintf(1, 'Opened SN ')
magic = fread(fid, 1, 'int32');
if (magic == hex2dec('1e3d4c51'))
	fprintf(1, 'binary')
	ndim = fread(fid, 1, 'int32');
else
	error('magic number not recognized')
end

dim = zeros(ndim);
fprintf(1, 'ndim = %d,  format = ', ndim);
size = 1;
for i=1:ndim,
  dim(i) = fread(fid, 1, 'int32');
  fprintf(1, '%d', dim(i))
  if  i < ndim
    fprintf(1, ' x ')
  end
  size = dim(i) * size;
end

fprintf(1, '\n')
while i<3
	if (magic ~= hex2dec('2e4d4154'))	% not for ascii format
		fread(fid, 1, 'int32');
	end
	i = i + 1;
end

if ndim == 1
	arr = zeros(dim(1),1);
elseif ndim == 2
	arr = zeros(dim(1), dim(2));
else
	arr = zeros(size,1);
	disp('I cannot define matrices with ndim>2, so I will just return it as a vector.')
end

if (magic == hex2dec('1e3d4c51'))
	if ndim == 2   % do this separately - reshaping afterwards takes memory
		for i=1:dim(1),
			arr(i,:) = fread(fid, dim(2), 'float')';
		end
	else
		arr = fread(fid, size, 'float')';
	end
end
fclose(fid);

mat = arr;




