function label_map = prep_label_check(CFG)

% prepare_label(CFG)
%
% Converts a genome annotation into a position-wise labeling used to
% train mTIM.
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009
%            Nico Goenitz, TU Berlin, 2013

CHECK = 1,
if CHECK,
  genome_info = init_genome(CFG.genome_info);
  for c=1:CFG.num_chr
    genome_info.flat_fnames{c} = [genome_info.basedir '/genome/' genome_info.contig_names{c} '.flat'];
    seq{c} = load_genomic(c,'+',1,CFG.chr_lens(c), genome_info);
  end
  cnt_good  = 0;
  cnt_total = 0;
end


% parse information from gene annotation
load(CFG.gene_fn, '-mat', 'genes');
LABELS = get_label_set_mTIM();
label_map = cell(1, CFG.num_chr);
for c=1:CFG.num_chr,
  label_map{c} = LABELS.intergenic*ones(1,CFG.chr_lens(c), 'int8');
end

tic
for g=1:length(genes),
  strand = genes(g).strand;
  chr = genes(g).chr_num;
  for t=1:length(genes(g).transcripts),
    exons = genes(g).exons{t};
    
    % this can happen in the Gencode annotations...
    if isequal(exons, [0 0]),
      continue
    end

    switch strand,
     case '+', 
      % convert half-open exon intervals into closed ones
      exons(:,2) = exons(:,2) - 1;
      for e=1:size(exons,1),
        idx = find(label_map{chr}(exons(e,1):exons(e,2)) ~= LABELS.intergenic ...
                   & label_map{chr}(exons(e,1):exons(e,2)) ~= LABELS.exon_W);

        idx = idx + exons(e,1) - 1;
        label_map{chr}(exons(e,1):exons(e,2)) = LABELS.exon_W;
        label_map{chr}(idx) = LABELS.ambiguous;
      end
      introns = [exons(1:end-1,2)+1, exons(2:end,1)-1];
      if CHECK, % check splice site consensus
        cnt_good = cnt_good + ...
            sum(seq{chr}(introns(:,1))   == 'g' & seq{chr}(introns(:,1)+1) == 't');
        cnt_good = cnt_good + ...
            sum(seq{chr}(introns(:,2)-1) == 'a' & seq{chr}(introns(:,2))   == 'g');
        cnt_total = cnt_total + 2*size(introns,1);
      end
      for i=1:size(introns,1),
        idx = find(label_map{chr}(introns(i,1):introns(i,2)) ~= LABELS.intergenic ...
                   & label_map{chr}(introns(i,1):introns(i,2)) ~= LABELS.intron_W);
        idx = idx + introns(i,1) - 1;
        label_map{chr}(introns(i,1):introns(i,2)) = LABELS.intron_W;
        label_map{chr}(idx) = LABELS.ambiguous;
      end
     case '-', 
      % convert half-open exon intervals into closed ones
      exons(:,1) = exons(:,1) + 1;
      for e=1:size(exons,1),
        idx = find(label_map{chr}(exons(e,1):exons(e,2)) ~= LABELS.intergenic ...
                   & label_map{chr}(exons(e,1):exons(e,2)) ~= LABELS.exon_C);
        
        idx = idx + exons(e,1) - 1;
        label_map{chr}(exons(e,1):exons(e,2)) = LABELS.exon_C;
        label_map{chr}(idx) = LABELS.ambiguous;
      end
      introns = [exons(1:end-1,2)+1, exons(2:end,1)-1];
      if CHECK, % check splice site consensus
        cnt_good = cnt_good + ...
            sum(seq{chr}(introns(:,1))   == 'c' & seq{chr}(introns(:,1)+1) == 't');
        cnt_good = cnt_good + ...
            sum(seq{chr}(introns(:,2)-1) == 'a' & seq{chr}(introns(:,2))   == 'c');
        cnt_total = cnt_total + 2*size(introns,1);
      end
      for i=1:size(introns,1),
        idx = find(label_map{chr}(introns(i,1):introns(i,2)) ~= LABELS.intergenic ...
                   & label_map{chr}(introns(i,1):introns(i,2)) ~= LABELS.intron_C);
        idx = idx + introns(i,1) - 1;
        label_map{chr}(introns(i,1):introns(i,2)) = LABELS.intron_C;
        label_map{chr}(idx) = LABELS.ambiguous;
      end
     otherwise error('unkown strand: %s', strand);
    end
  end
  if mod(g,100)==0,
    fprintf('  processed %i genes (%2.1f%%, %.1f sec)\r', g, 100*g/length(genes), toc);
  end
end
fprintf('  processed %i genes (%2.1f%%, %.1f sec)\n', g, 100*g/length(genes), toc);

if CHECK,
  fprintf('  %2.2f%% of all splice sites had consensus GT-AG\n', ...
          100*cnt_good/cnt_total );
end


% add overlapping genes ambiguous blocks
cnt = 0;
for c=1:CFG.num_chr,
  inds = find([genes.chr_num] == c);

  lens = [genes.stop]-[genes.start] + 1;
  offs = [genes.start];
  ends = [genes.stop];

  
  for i=1:length(inds),
      % rm i from inds
      ind = inds(i);
      ninds = setdiff(inds,ind);

      % left side (or completely) overlapping 
      ainds = find(offs(ninds) <= offs(ind) ...
                    & ends(ninds) >= offs(ind) ...
                    & ends(ninds) <= ends(ind));
      if ~isempty(ainds),
      
        for j=1:length(ainds),
            idx = ninds(ainds(j));
            amb_block = offs(ind):ends(idx);
            cnt = cnt + sum(label_map{c}(amb_block) ~= LABELS.ambiguous);
            label_map{c}(amb_block) = LABELS.ambiguous;
        end
      end

      % right side overlapping
      ainds = find(offs(ninds) <= ends(ind) ...
                    & offs(ninds) > offs(ind) ...
                    & ends(ninds) > ends(ind));
      if ~isempty(ainds),
      
        for j=1:length(ainds),
            idx = ninds(ainds(j));
            amb_block = offs(idx):ends(ind);
            cnt = cnt + sum(label_map{c}(amb_block) ~= LABELS.ambiguous);
            label_map{c}(amb_block) = LABELS.ambiguous;
        end
      end

      % both sides overlapping
      ainds = find(offs(ninds) < offs(ind) ...
                    & ends(ninds) > ends(ind));
      if ~isempty(ainds),
      
        for j=1:length(ainds),
            idx = ninds(ainds(j));
            amb_block = offs(ind):ends(ind);
            cnt = cnt + sum(label_map{c}(amb_block) ~= LABELS.ambiguous);
            label_map{c}(amb_block) = LABELS.ambiguous;
        end
      end

      % completely inside
      ainds = find(offs(ninds) > offs(ind) ...
                    & ends(ninds) < ends(ind));
      if ~isempty(ainds),
      
        for j=1:length(ainds),
            idx = ninds(ainds(j));
            amb_block = offs(idx):ends(idx);
            cnt = cnt + sum(label_map{c}(amb_block) ~= LABELS.ambiguous);
            label_map{c}(amb_block) = LABELS.ambiguous;
        end
      end
  end

end
fprintf('Labeled %i additional ambiguous positions.\n',cnt);


fn = fieldnames(LABELS);
for f=1:length(fn),
  cnt = 0;
  for c=1:CFG.num_chr,
    cnt = cnt + sum(label_map{c} == getfield(LABELS, fn{f}));
  end
  fprintf('Labeled %12i genomic positions as %10s (%2.1f%%)\n', ...
          cnt, fn{f}, 100*cnt/sum(CFG.chr_lens));
end 


%save(CFG.label_fn, '-v7.3', 'CFG', 'LABELS', 'label_map');


% eof
