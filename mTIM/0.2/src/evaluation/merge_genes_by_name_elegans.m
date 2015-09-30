function genes = merge_genes_by_name_elegans(genes)
% genes = merge_genes_by_name_elegans(genes)
% 
% Merges transcripts with the same gene name into a single gene structure.
%
% written by Gunnar Raetsch & Georg Zeller, MPI Tuebingen, 2007-2008

fprintf('Merging by name, starting with %i genes\n', length(genes));

gnames = {genes.name};
% this is C. elegans-specific!
for g = 1:length(gnames),
  gnames{g} = gnames{g}(1:find(gnames{g}=='.',1,'last')-1);
end

[tmp, tmp, idx] = unique(gnames);
[list, idx2] = sort(idx);
take_map = ones(1,length(genes));
for i=1:length(idx)-1,
  j = 0;
  merge_occurred = 1;
  while merge_occurred,
    merge_occurred = 0;
    j = j + 1;
    if i+j<=length(list) && list(i)==list(i+j), 
      % only merge if genes are on the same strand and overlapping
      if genes(idx2(i)).strand == genes(idx2(i+j)).strand ...
            && genes(idx2(i)).start <= genes(idx2(i+j)).stop ...
            && genes(idx2(i)).stop >= genes(idx2(i+j)).start,
        merge_occurred = 1;
        assert(isequal(genes(idx2(i)).chr, genes(idx2(i+j)).chr));
        
        genes(idx2(i)).transcripts = {genes(idx2(i)).transcripts{:} ...
                            genes(idx2(i+j)).transcripts{:}};
        genes(idx2(i)).transcript_valid = [genes(idx2(i)).transcript_valid ...
                            genes(idx2(i+j)).transcript_valid];
        
        genes(idx2(i)).exons = {genes(idx2(i)).exons{:} ...
                          genes(idx2(i+j)).exons{:}};
        
        assert(length(genes(idx2(i)).transcripts) ...
               == length(genes(idx2(i)).exons));
        
        % cds_exons and UTRs
        if isfield(genes, 'cds_exons'),
          genes(idx2(i)).cds_exons = {genes(idx2(i)).cds_exons{:} ...
                              genes(idx2(i+j)).cds_exons{:}};
        end
        if isfield(genes, 'utr5_exons'),
          genes(idx2(i)).utr5_exons = {genes(idx2(i)).utr5_exons{:} ...
                            genes(idx2(i+j)).utr5_exons{:}};
        end
        if isfield(genes, 'utr3_exons'),
          genes(idx2(i)).utr3_exons = {genes(idx2(i)).utr3_exons{:} ...
                              genes(idx2(i+j)).utr3_exons{:}};
        end
        
        % gene properties
        if isfield(genes, 'is_alt'),
          genes(idx2(i)).is_alt = any([genes(idx2(i)).is_alt, ...
                              genes(idx2(i+j)).is_alt]);
        end
        if isfield(genes, 'is_alt_spliced'),
          genes(idx2(i)).is_alt_spliced = any([genes(idx2(i)).is_alt_spliced, ...
                            genes(idx2(i+j)).is_alt_spliced]);
        end
        if isfield(genes, 'is_valid'),
          genes(idx2(i)).is_valid = all([genes(idx2(i)).is_valid, ...
                              genes(idx2(i+j)).is_valid]);
        end
        if isfield(genes, 'transcript_complete'),
          genes(idx2(i)).transcript_complete = [genes(idx2(i)).transcript_complete, ...
                              genes(idx2(i+j)).transcript_complete];
        end
        if isfield(genes, 'transcript_coding'),
          genes(idx2(i)).transcript_coding = [genes(idx2(i)).transcript_coding, ...
                              genes(idx2(i+j)).transcript_coding];
        end
        
        take_map(idx2(i+j)) = 0;
      end
    end
  end
  start = inf; 
  stop  = 0;
  for j=1:length(genes(i).transcripts),
    if start>min(genes(i).exons{j}(:)), 
      start = min(genes(i).exons{j}(:));
    end
    if stop<max(genes(i).exons{j}(:)),
      stop = max(genes(i).exons{j}(:)); 
    end
  end
  assert(~isinf(start) & stop~=0);
  genes(i).start = start;
  genes(i).stop  = stop;
end
idx = find(take_map);
genes = genes(idx);
fprintf('%i genes remaining after merging by name\n', length(genes));
