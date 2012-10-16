function genes = merge_transcripts_by_colocation(genes, fields, merge_by_cds, verb)
% MERGE_TRANSCRIPTS_BY_COLOCATION   Merges transcripts from overlapping gene loci.
%
%   genes = merge_transcripts_by_colocation(genes, fields, merge_by_cds, verb)
%
%   -- input --
%   genes:        struct defining genes with start, stops, exons etc.
%   fields:       fields of genes struct that should be considered
%                 when merging
%   merge_by_cds: merge transcripts according to their CDS exons
%   verb:         verbosity
%
%   -- output --
%   genes:        struct defining genes with start, stops, exons etc.,
%                 where overlapping transcripts are merged
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2010-2011 Jonas Behr, Gunnar Raetsch, Regina Bohnert
%   Copyright (C) 2010-2011 Max Planck Society
%


if nargin==1 || isempty(fields)
  fields = {'transcripts', 'exons', 'cds_exons', 'utr3_exons', 'utr5_exons', 'tss', 'tis', 'cdsStop', 'cleave'};
end

if nargin<3
  merge_by_cds = 0;
end

[genes.parent_gene] = deal(0);
if ~isfield(genes, 'strands')
  [genes.strands] = deal('');
end
for g = 1:length(genes),
  genes(g).parent_gene = g*ones(1, length(genes(g).transcripts));
  genes(g).strands = repmat(genes(g).strand, 1, length(genes(g).transcripts));
end
fields = {fields{:}, 'parent_gene', 'strands'};

idx = 1;
while ~isempty(idx)
  nofgenes = length(genes);
  genes = find_start_stop(genes, merge_by_cds);
  [genes.remove] = deal(0);
  idx = find_overlapping_regions(genes, [], 0);
  if exist('issorted')
    assert(issorted(idx(:, 1)));
  end
  for row = 1:size(idx, 1),
    % avoid duplicating transcripts
    if genes(idx(row, 1)).remove==1 || genes(idx(row, 2)).remove==1; continue; end
    genes(end+1) = merge_genes(genes(idx(row, 1)), genes(idx(row,2)), fields);
    genes(idx(row, 1)).remove = 1;
    genes(idx(row, 2)).remove = 1;
  end
  rm_idx = find([genes.remove]);
  if nargin<4 || verb
    fprintf(1, 'reducing number of genes from %i to %i\n', nofgenes, nofgenes-length(rm_idx)/2);
  end
  genes(rm_idx) = [];
  if isempty(rm_idx)
    break; 
  end
end
genes = find_start_stop(genes, 0);
genes = rmfield(genes, 'remove');
% update is_alt field
for g = 1:length(genes),
  if length(genes(g).transcripts)>1
    genes(g).is_alt = 1;
  else
    genes(g).is_alt = 0;
  end
end
return;


% adapts start and stop coordinates of genes
function genes = find_start_stop(genes, merge_by_cds)
for g = 1:length(genes),
  start = inf; stop = 0;
  if merge_by_cds && isfield(genes, 'cds_exons') && ~isempty(genes(g).cds_exons)
    for t = 1:length(genes(g).cds_exons),
      start = min(start, min(genes(g).cds_exons{t}(:)));
      stop = max(stop, max(genes(g).cds_exons{t}(:)));
    end
  else
    for t = 1:length(genes(g).exons),
      start = min(start, min(genes(g).exons{t}(:)));
      stop = max(stop, max(genes(g).exons{t}(:)));
    end
  end
  assert(~isinf(start) & stop~=0);
  genes(g).start = start;
  genes(g).stop = stop;
end
return;


% merges two genes and their corresponding fields
function gene1 = merge_genes(gene1, gene2, fields)
for f = 1:length(fields),
  if isfield(gene1, fields{f})
    if iscell(gene1.(fields{f})) && iscell(gene2.(fields{f}))
      gene1.(fields{f}) = {gene1.(fields{f}){:} gene2.(fields{f}){:}};
    elseif iscell(gene1.(fields{f}))
      gene1.(fields{f}) = {gene1.(fields{f}){:} {}};
    elseif iscell(gene1.(fields{f}))
      gene1.(fields{f}) = {{} gene2.(fields{f}){:}};
    else
      gene1.(fields{f}) = [gene1.(fields{f}) gene2.(fields{f})];
    end
  end
end 
return;

