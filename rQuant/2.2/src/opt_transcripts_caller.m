function genes = opt_transcripts_caller(PAR)
% OPT_TRANSCRIPTS_CALLER   Prepares data for transcript weight optimisation problem.
%
%   genes = opt_transcripts_caller(PAR)
%
%   -- input --
%   PAR contains
%     CFG:             configuration struct
%     genes:           struct defining genes with start, stops, exons etc.
%     profile_weights: weights of profile functions
%
%   -- output --
%   genes:             struct with additional fields of e.g. estimated transcript weights
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


CFG = PAR.CFG;
genes = PAR.genes;
profile_weights = PAR.profile_weights;
clear PAR;

%%%% paths
addpath(CFG.paths);

[genes.mean_ec] = deal([]);
[genes.coverage] = deal([]);
[genes.introns] = deal([]);
[genes.transcript_weights] = deal([]);
[genes.obj] = deal([]);

chr_num = unique([genes.chr_num]);
for c = chr_num,
  chr_idx = find([genes.chr_num]==c);
  if CFG.VERBOSE>0, fprintf(1, '\nprocessing %i genes for contig %i (%s)\n', length(chr_idx), genes(chr_idx(1)).chr_num, genes(chr_idx(1)).chr); end
  for g = chr_idx,
    gene = genes(g);
    if CFG.VERBOSE>0, fprintf(1, '\ngene %i: %i isoform(s) with %i exonic positions\n', g, length(gene.transcripts), gene.exonic_len); end
    if 0
      if CFG.VERBOSE>0, fprintf(1, 'gene too complex %i\n', g); end
      genes(g).transcript_weights(1:length(genes(g).transcripts)) = nan;
      genes(g).obj = nan;
      continue;
    end
    %%%%% load exon coverage and introns for gene %%%%%
    try
      if CFG.VERBOSE>1, fprintf(1, 'Loading reads...\n'); tic; end
      if CFG.paired
        [coverage, reads_ok, introns, read_starts, paired_reads] = get_coverage_per_read(CFG, gene);
        clear read_starts;
      else
        [coverage, reads_ok, introns] = get_coverage_per_read(CFG, gene);
      end
      if CFG.VERBOSE>1, fprintf(1, 'Took %.1fs.\n', toc); end
    catch
      reads_ok = 0;
    end
    if ~reads_ok,
      if CFG.VERBOSE>0, fprintf(1, 'coverage could not be loaded for gene %i\n', g); end
      genes(g).transcript_weights(1:length(genes(g).transcripts)) = nan;
      genes(g).obj = nan;
      continue;
    end
    genes(g).mean_ec = full(mean(sum(coverage,2)));
    genes(g).coverage = sum(coverage,2);
    genes(g).introns = introns;
    if isempty(coverage) || genes(g).mean_ec==0
      if CFG.VERBOSE>0, fprintf(1, 'no coverage for gene %i\n', g); end
      genes(g).transcript_weights(1:length(genes(g).transcripts)) = 0;
      genes(g).obj = 0;
      continue;
    end
    %%%%% prepare exon mask %%%%%
    exon_mask = zeros(gene.exonic_len, length(gene.transcripts));
    for t = 1:length(gene.transcripts),
      if isfield(gene, 'strands') && ~isempty(gene.strands),
        strand_str = gene.strands(t);
      else
        strand_str = gene.strand;
      end
      % fill exon mask
      if strand_str=='+'
        rev_idx = 1:size(profile_weights,1);
      else
        rev_idx = size(profile_weights,1):-1:1;
      end
      % evaluate profile PLiFs for each exonic position
      [feat feat_val feat_val_next] = gen_exon_features(CFG, gene, t);
      fidx = find(feat_val>0); 
      feat = feat(fidx,:); feat_val = feat_val(fidx,:); feat_val_next = feat_val_next(fidx,:);
      exon_mask(fidx,t) = gen_exon_mask(profile_weights(rev_idx,:), gene.transcript_len_bin(t), feat, feat_val, feat_val_next, [1:length(fidx)]', ones(length(fidx),1));
    end
    
    %%%%% prepare intron mask %%%%%
    if ~isempty(introns)
      [intron_mask intron_count] = get_intron_data(gene, CFG, introns, g);
    else
      intron_mask = zeros(0,length(gene.transcripts));
      intron_count = [];
    end
      
    %%%%% prepare repeat mask %%%%%
    repeat_mask = false(gene.exonic_len, 1); 
    fname = sprintf('%s/%s_repeat', CFG.repeats_fn, gene.chr);
    if exist(sprintf('%s.pos', fname), 'file')
      [map.pos map.repeats] = interval_query(fname, {'repeats'}, [gene.start;gene.stop]);
      if ~isempty(map.pos)
        [tmp idx1 idx2] = intersect(map.pos, gene.eidx);
        assert(length(idx2)<=length(map.pos));
        repeat_mask(idx2) = true;
      end
    end
    
    %%%%% segment or position-based %%%%%
    if strcmp(CFG.method, 'seg') % segment-wise
      % code only works for CFG.paired=0
      assert(size(coverage,1)==gene.exonic_len);
      segments = gen_segments(gene);
      assert(sum(segments(:,2)-segments(:,1)+1)==gene.exonic_len);
      exon_mask_pos = exon_mask;
      coverage_pos = coverage;
      exon_mask = zeros(size(segments,1),length(gene.transcripts));
      coverage = zeros(size(segments,1),1);
      for s = 1:size(segments,1),
        sidx = segments(s,1):segments(s,2);
        exon_mask(s,:) = median(exon_mask_pos(sidx,:),1);
        coverage(s) = round(median(coverage_pos(sidx)));
      end
      clear exon_mask_pos coverage_pos;
    elseif strcmp(CFG.method, 'pos') % position-wise
      % only consider regions in at least one exon
      if size(exon_mask,1)<2*(CFG.max_side_len-1)
        if ~all(any(exon_mask, 2)')
          [mval midx] = max(gene.transcript_len_bin);
          if CFG.VERBOSE>0, fprintf('gene inconsistent with profile (%i: %i)\n', mval, gene.transcript_length(midx)); end
        end
        %assert(all(any(exon_mask, 2)'));
      else
        %assert(all(any(exon_mask([1:CFG.max_side_len-1, end-CFG.max_side_len+2:end],:),2)'));
      end  
      mask = (any(exon_mask, 2) & ~repeat_mask)';
      exon_mask = exon_mask(mask, :);
      coverage = coverage(mask, :);
    else
      error('unknown method %s', CFG.method);
    end
    if isempty(exon_mask)
      if CFG.VERBOSE>0, fprintf(1, 'no positions left for gene %i\n', g); end
      genes(g).transcript_weights(1:length(genes(g).transcripts)) = nan;
      genes(g).obj = nan;
      continue;
    end
    
    %%%%% generate paired-end data %%%%%
    if CFG.paired
      [paired_exp, paired_obs] = get_paired_data(CFG, gene, paired_reads);
    else
      paired_exp = []; paired_obs = [];
    end
    
    %%%%% transcript weight optimisation %%%%%
    C_w = gene.transcript_length';
    [weights, obj] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, intron_mask, paired_exp, paired_obs, C_w);
    genes(g).transcript_weights = weights;
    genes(g).obj = obj;
  end
end
