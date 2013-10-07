function savesn(MATRIX, FILENAME, FORMAT)
% SAVESN save a matrix in SN file format
%	SAVESN (MATRIX, FILENAME, FORMAT) saves MATRIX as an SN matrix in a
%       file with name FILENAME. By default, '.mat' is added to the FILENAME.
%	FORMAT can be 'binary' (default), 'double', 'ascii', or 'packed'.

if nargin == 3
	if strcmp(FORMAT, 'ascii')
		magic = hex2dec('2e4d4154');
	elseif strcmp(FORMAT, 'binary')
		magic = hex2dec('1e3d4c51');
	elseif strcmp(FORMAT, 'double')
		magic = hex2dec('1e3d4c53');
	elseif strcmp(FORMAT, 'packed')
		magic = hex2dec('1e3d4c52');
	else
		error('format not recognized')
	end
else
	magic = hex2dec('1e3d4c51');
end

if findstr(FILENAME, '.')    % check for extension .mat, .pat, .nam or other
	[fid, message] = fopen(FILENAME, 'w');
else
	[fid, message] = fopen([FILENAME '.mat'], 'w');
end

% as vectors in matlab are nx1 or 1xn matrices, we save them with ndim=2, too.
% Thus, information is preserved about whether it's a row or a column vector.
ndim = 2;
dim = size(MATRIX);

fprintf(1, '%s', message)
fprintf(1, 'Created SN ')
fwrite(fid, magic, 'int32');
if (magic == hex2dec('2e4d4154'))
	fprintf(1, 'ascii')
	fprintf(fid, ' %d', ndim);
elseif (magic == hex2dec('1e3d4c51'))
	fprintf(1, 'binary')
	fwrite(fid, ndim, 'int32');
elseif (magic == hex2dec('1e3d4c53'))
	fprintf(1, 'double')
	fwrite(fid, ndim, 'int32');
elseif (magic == hex2dec('1e3d4c52'))
	fprintf(1, 'packed')
	fwrite(fid, ndim, 'int32');
else
	error('magic number not recognized')
end
fprintf(1, ' matrix file, ')

fprintf(1, 'ndim = %d,  format = ', ndim);
size = 1;
for i=1:ndim,
	if (magic == hex2dec('2e4d4154'))	% ascii format
		fprintf(fid, ' %d', dim(i));
	else
		fwrite(fid, dim(i), 'int32');
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
		fwrite(fid, 1, 'int32');	% SN1 compatibility
	end
	i = i + 1;
end
if (magic == hex2dec('2e4d4154'))	% ascii format
		fprintf(fid, '\n');
end

fprintf(1, 'writing data ');
if (magic == hex2dec('1e3d4c51'))	%binary
	dots_ctr = 0;
	dot_every = dim(1)/10;
	for i=1:dim(1),
		dots_ctr = dots_ctr + 1;
		if dots_ctr > dot_every
			fprintf(1, '.');
			dots_ctr = 0;
		end
		fwrite(fid, MATRIX(i,:)', 'float');
	end
	fprintf(1, '\n');
end

fprintf(1, 'writing data ');
if (magic == hex2dec('1e3d4c53'))	%binary double
	dots_ctr = 0;
	dot_every = dim(1)/10;
	for i=1:dim(1),
		dots_ctr = dots_ctr + 1;
		if dots_ctr > dot_every
			fprintf(1, '.');
			dots_ctr = 0;
		end
		fwrite(fid, MATRIX(i,:)', 'double');
	end
	fprintf(1, '\n');
end

if (magic == hex2dec('1e3d4c52'))	%packed
	for i=1:dim(1),
		for j=1:dim(2),
			if MATRIX(i,j) >= 8
				b = hex2dec('7f');
				fprintf(1,'clipped value %f to 8 in packed saving\n', MATRIX(i,j));
			elseif MATRIX(i,j) <= -8
				b = hex2dec('80');
				fprintf(1,'clipped value %f to -8 in packed saving\n', MATRIX(i,j));
			else %if MATRIX(i,j) >= 0;
				b = MATRIX(i,j) * 16;
%			else
%				b = (b - 256) / 16;
			end
			fwrite(fid, b, 'unsigned char');
		end
	end
end
if (magic == hex2dec('2e4d4154'))	%ascii
	for i=1:dim(1),
		 fprintf(fid, '%f\n', MATRIX(i,:)');
	end
end
fclose(fid);




