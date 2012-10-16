function [intron_mask intron_count] = get_intron_data(gene, CFG, introns, gene_idx)
% GET_INTRON_DATA   Prepares intron data for gene.
%
%   [intron_mask intron_count] = get_intron_data(gene, CFG, introns)
%
%   -- input --
%   gene:         struct defining a gene with start, stops, exons etc.
%   CFG:          configuration struct
%   introns:      nx4 list of introns (intron start, intron stop, confirmation, strand)
%   gene_idx:     index of gene
%
%   -- output --
%   intron_mask:  IxT matrix defining whether an intron belongs to a particular transcript
%   intron_count: vector of observed intron confirmation
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2011 Max Planck Society
%


% plus strand
fidx = find(introns(:,4)==0);
intron_starts{1} = introns(fidx,1); 
intron_stops{1} = introns(fidx,2);
conf{1} = introns(fidx,3);
% minus strand
fidx = find(introns(:,4)==1);
intron_starts{2} = introns(fidx,1); 
intron_stops{2} = introns(fidx,2);
conf{2} = introns(fidx,3);

%%%%% prepare intron list %%%%%
intron_list = zeros(0,4); 
num_transcripts(1) = length(gene.transcripts);
if CFG.both_strands && isfield(gene, 'strands') && ~isempty(gene.strands)
  strand_str = unique(gene.strands);
  num_transcripts(2) = sum(gene.strands=='+');
  num_transcripts(3) = sum(gene.strands=='-');
else
  strand_str = gene.strand;
  if gene.strand=='+',
    num_transcripts(2) = length(gene.transcripts);
    num_transcripts(3) = 0;
  else
    num_transcripts(2) = 0;
    num_transcripts(3) = length(gene.transcripts);
  end
end
assert(length(strand_str)==1|length(strand_str)==2);
for s = 1:length(strand_str),
  strand = (strand_str(s)=='-') + 1;
  if ~isempty(intron_starts{strand})
    idx = find(intron_starts{strand}>=gene.start & intron_stops{strand}<=gene.stop);
    if ~isempty(idx)
      intron_list = [intron_list; ...
                     double([intron_starts{strand}(idx), intron_stops{strand}(idx), ...
                     conf{strand}(idx), ((strand_str(s)=='-')+1)*ones(length(idx),1)])];
    end
  end
end
if CFG.VERBOSE>=2,
  fprintf(1, 'found %i transcripts, %i on + strand, %i on - strand\n', num_transcripts(1), num_transcripts(2), num_transcripts(3));
  if exist('read_starts', 'var')
    fprintf(1, 'found %i reads\n', sum(read_starts));
  end
  fprintf(1, 'found %i introns, %i on + strand, %i on - strand (%i ignored)\n', size(intron_list,1), sum(intron_list(:,4)==1), sum(intron_list(:,4)==2), size(introns,1)-size(intron_list,1));
end

intron_mask = zeros(size(intron_list,1), length(gene.transcripts));
for t = 1:length(gene.transcripts),
  % fill intron mask
  exons = gene.exons{t};
  num_found = 0;
  for e = 1:size(exons,1)-1,
    if ~isempty(intron_list)
      idx = find(exons(e,2)+1==intron_list(:,1) & exons(e+1,1)-1==intron_list(:,2));
      if CFG.both_strands && isfield(gene, 'strands') && ~all(intron_list(idx,4)==(gene.strands(t)=='-')+1)
        idx = [];
      end
      if CFG.both_strands && isfield(gene, 'strands') && ~isempty(idx)
        assert(all(intron_list(idx,4)==(gene.strands(t)=='-')+1));
      end
    end
    if ~isempty(intron_mask)
      if ~isempty(idx),
        num_found = num_found + 1;
      end
      intron_mask(idx,t) = 1;
    end
  end
  if num_found==0,
    if CFG.VERBOSE>0
      fprintf(1, 'introns not found in gene %i, transcript %i \n', gene_idx, t);
    end
  else
    if CFG.VERBOSE>=2
      fprintf(1, 'found %i matching introns, gene %i, transcript %i\n', num_found, gene_idx, t)
    end
  end
end

    
if ~isempty(intron_mask)
  intron_count = intron_list(:,3);
else
  intron_count = [];
end