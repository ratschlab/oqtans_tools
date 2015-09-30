function match = match_transcripts(gene, ovl_genes)

% match will be 1 if any transcript among those in ovl_genes matches a
% transcript of gene
assert(length(gene) == 1);

gene_transcripts = get_transcript_introns(gene);
ovl_transcripts = get_transcript_introns(ovl_genes);

if length(gene_transcripts)*length(ovl_transcripts) > 100,
  warning('many transcripts to compare!')
%  keyboard
end

match = 0;
for g=1:length(gene_transcripts),
  for o=1:length(ovl_transcripts),
    if isequal(gene_transcripts{g}, ovl_transcripts{o}),
      match = 1;
      break
    end
  end
end

%fprintf('gene transcripts\n'); 
%gene_transcripts{:}
%fprintf('overlapping transcripts\n'); 
%ovl_transcripts{:}
%match
%keyboard
