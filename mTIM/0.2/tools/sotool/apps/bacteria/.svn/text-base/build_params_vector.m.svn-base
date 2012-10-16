function [mask] = build_params_vector(params)
% -- description --
% Builds a matrix containing all unique combinations of 
% states that the params-cell-vectors defined.
%
% -- in --
%   params  : cell array of vectors (=parameters) containing all
%             possible states of the current parameter
%
% -- out --
%   mask : Returns a matrix with #parameters columns
%
% -- usage --
%   params = {[0,1],[0,1],[0,1]};
%   flattenParams(params) = [0,0,0;
%                            0,0,1;
%                            0,1,0;
%                               ...
%
% author: nico goernitz

% number of parameters
paramCount = length(params);

% start with a non-empty mask
mask = params{1}';

% for every parameter
for i=2:paramCount
    if (~isempty(params{i})),
      % append all states of current parameter
      % for every line in 'mask'
      mask = append(mask, params{i});
    end
end


function m = append(mask,p)
% Append column p for every line in m.
% - p has to be a column vector
% - mask is a matrix
m = [];
for i=1:size(mask,1)
    for j=1:length(p)
        value = [mask(i,:), p(j)];
        m = [m; value];
    end
end