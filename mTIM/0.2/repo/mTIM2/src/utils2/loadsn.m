function [mat] = loadsn(FILENAME)
% LOADSN load a matrix in SN file format
%	LOADSN (FILENAME) loads an SN matrix (binary or ascii) from
%	the file with name FILENAME

if findstr(FILENAME, '.')    % check for extension .mat, .pat, .nam or other
	[fid, message] = fopen(FILENAME);
else
	[fid, message] = fopen([FILENAME '.mat']);
end

fprintf(1, '%s', message)
fprintf(1, 'Opened SN ')
magic = fread(fid, 1, 'int32');
if (magic == hex2dec('2e4d4154'))
	fprintf(1, 'ascii')
	ndim = fscanf(fid, ' %d', 1);
elseif (magic == hex2dec('1e3d4c51'))
	fprintf(1, 'binary')
	ndim = fread(fid, 1, 'int32');
elseif (magic == hex2dec('1e3d4c53'))
	fprintf(1, 'double')
	ndim = fread(fid, 1, 'int32');
elseif (magic == hex2dec('1e3d4c52'))
	fprintf(1, 'packed')
	ndim = fread(fid, 1, 'int32');
else
	error('magic number not recognized')
end
fprintf(1, ' matrix file, ')

dim = zeros(ndim);
fprintf(1, 'ndim = %d,  format = ', ndim);
size = 1;
for i=1:ndim,
	if (magic == hex2dec('2e4d4154'))	% ascii format
		dim(i) = fscanf(fid, ' %d', 1);
	else
		dim(i) = fread(fid, 1, 'int32');
	end
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
if (magic == hex2dec('1e3d4c53'))
	if ndim == 2   % do this separately - reshaping afterwards takes memory
		for i=1:dim(1),
			arr(i,:) = fread(fid, dim(2), 'double')';
		end
	else
		arr = fread(fid, size, 'float')';
	end
end
if (magic == hex2dec('1e3d4c52'))
	if ndim == 2
		for i=1:dim(1),
			arr(i,:) = fread(fid, dim(2), 'unsigned char')';
			for j=1:dim(2),
				if arr(i,j) == hex2dec('7f')
					arr(i,j) = 8;
				elseif arr(i,j) == hex2dec('80')
					arr(i,j) = -8;
				elseif arr(i,j) < hex2dec('7f')
					arr(i,j) = arr(i,j)/16;
				elseif arr(i,j) > hex2dec('80')
					arr(i,j) = (arr(i,j) - 256) / 16;
				end
			end
		end
	else
		arr = fread(fid, size, 'unsigned char')';
        	for i=1:size,
			if arr(i) == hex2dec('7f')
				arr(i) = 8;
			elseif arr(i) == hex2dec('80')
				arr(i) = -8;
			elseif arr(i) < hex2dec('7f')
				arr(i) = arr(i)/16;
			elseif arr(i) > hex2dec('80')
				arr(i) = (arr(i) - 256) / 16;
			end
		end
	end
end
if (magic == hex2dec('2e4d4154'))
	if ndim == 2
		for i=1:dim(1),
			arr(i,:) = fscanf(fid, '%f\n', dim(2))';
		end
	else
		arr = fscanf(fid, '%f\n', size)';
	end
end
fclose(fid);

mat = arr;
