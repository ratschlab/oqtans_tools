function txt = print_structure(IN,space)
% This little piece of code will take
% as input any structure and plot its
% fields and values in a pretty format.
% Should any of the subfields be a
% structure as well then the function will
% be called recursivly.
% 
%
%

tab = '';
if exist('space'),
    tab = space;
end


txt = '';
names = fieldnames(IN);
names = sort(names);

for i=1:length(names),
    cont = getfield(IN,names{i});

   if isstruct(cont),
       txt = [txt sprintf('%s%s:\n',tab,names{i})];
       txt = [txt print_structure(cont,[tab '   '])]; 
   else
       val = get_value(cont);
       txt = [txt sprintf('%s%s: %s\n',tab,names{i},val)];
   end

end



function out = get_value(val)

out = '(not defined)';

if isa(val,'function_handle'),
    out = sprintf('%s (method)',func2str(val));
    return;
end

if isempty(val),
    out = '(emtpy)';
    return;
end

if iscell(val),
    out = '(cell)';
    return;
end

if isscalar(val),
    if isinteger(val),
        out = sprintf('%i',val);
    else
        out = sprintf('%1.4f',val);
    end
    return;
end

if ischar(val),
    out = sprintf('%s',val);
    return;
end


