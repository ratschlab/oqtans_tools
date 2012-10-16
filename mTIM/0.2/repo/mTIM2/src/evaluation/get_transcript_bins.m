function levels = get_transcript_bins(genes, num_bins)
% Sort genes according to their respective amount of transcripts.
% If 'num_bins' is present discretize them into 'num_bins' bins else
% use for each number of transcript a single bin resulting in highly
% unbalanced and even empty bins.
%
% IN
%   genes     : [1..n] gene structures
%   num_bins  : (optional) number of bins. if present the genes are
%               descritized according to their respective amount of
%               transcripts.
%
%
% OUT
%  levels     : [1..num_bins] cell array. each entry contains a list
%               with the associated gene indices.
%
% written by nico goernitz


% sort all genes according to number of transcripts
% each transcript is a bin
levels = {};
transcripts = zeros(1,length(genes));

for i=1:length(genes),
    % number of transcripts
    t = length(genes(i).transcripts);
    transcripts(i) = t;
    if length(levels)<t, levels{t} = [i];
    else
        if isempty(levels{t}), levels{t} = [i];
        else levels{t} = [levels{t}, i]; end
    end
end

if (exist('num_bins')),

  % TODO
  warning('NOT implemented yet.');

end
