function [pair_mat_exp, pair_mat_obs, segments, paired_reads_exp] = get_paired_data(CFG, gene, paired_reads_obs)
% GET_PAIRED_DATA   Prepares paired data for gene.
%
%   [pair_mat_exp, pair_mat_obs, segments, paired_reads_exp] = get_paired_data(CFG, gene, paired_reads_obs)
%
%   -- input --
%   CFG:              configuration struct
%   gene:             struct defining a gene with start, stops, exons etc.
%   paired_reads_obs: struct including starts, stops and mates for observed paired-end reads
%
%   -- output --
%   pair_mat_exp:     matrix (#segments x #segments x T) of connectivity by expected paired-end reads
%   pair_mat_obs:     matrix (#segments x #segments) of connectivity by observed paired-end reads
%   segments:         distinguishable segments
%   paired_reads_exp: struct including starts, stops and mates for expected paired-end reads
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2011 Regina Bohnert
%   Copyright (C) 2011 Max Planck Society
%


var_ins_size = true;
if isempty(CFG.ins_sizes)
  CFG.ins_sizes = 100;
end


%%% get distinguishable segments from splicegraph
segments = define_segments(gene.splicegraph{1}, gene.splicegraph{2});
len_segments = segments(:,2)-segments(:,1)+1;
pair_mat_len = zeros(size(segments,1));
for s1 = 1:size(segments,1),
 for s2 = 1:size(segments,1),
   if s1~=s2
     pair_mat_len(s1,s2) = len_segments(s1)+len_segments(s2);
   else
     pair_mat_len(s1,s2) = len_segments(s1);
   end
  end
end

%%% generate paired-end reads from annotated transcripts
pair_mat_exp = zeros(size(segments,1), size(segments,1), length(gene.transcripts));
for t = 1:length(gene.transcripts),
  tidx = [];
  for e = 1:size(gene.exons{t},1),
    tidx = [tidx, gene.exons{t}(e,1):gene.exons{t}(e,2)];
    tidx = unique(tidx);
  end
  if length(tidx) < 2*CFG.read_len+min(CFG.ins_sizes)
    continue;
  end
  if var_ins_size
    % insert size drawn from empirical distribution
    max_iter = min(5, length(CFG.ins_sizes));
    num_reads = 0; num_reads_exp = 2*max_iter*(length(tidx)-2*CFG.read_len+1);
    paired_reads_exp.starts = nan(1, num_reads_exp);
    paired_reads_exp.stops = nan(1, num_reads_exp);
    paired_reads_exp.mates = nan(2, num_reads_exp);
    for n = 1:length(tidx)-2*CFG.read_len+1,
      ridx = randperm(length(CFG.ins_sizes));
      ridx = ridx(1:max_iter);
      cand1_start = tidx(n);
      cand1_stop = tidx(n+CFG.read_len-1);
      for iter = 1:max_iter,
        ins_size = CFG.ins_sizes(ridx(iter)); % sampled insert size from observed distribution
        if n+2*CFG.read_len+ins_size>length(tidx)
          continue;
        end
        cand2_start = tidx(n+CFG.read_len+ins_size);
        cand2_stop = tidx(n+2*CFG.read_len+ins_size-1);
        paired_reads_exp.starts(num_reads+[1:2]) = [cand1_start, cand2_start];
        paired_reads_exp.stops(num_reads+[1:2]) = [cand1_stop, cand2_stop];
        paired_reads_exp.mates(:,num_reads/2+1) = [num_reads+1; num_reads+2];
        num_reads = num_reads + 2;
      end
    end
    assert(num_reads<=num_reads_exp);
    if num_reads>0
      paired_reads_exp.starts = paired_reads_exp.starts(1:num_reads);
      paired_reads_exp.stops = paired_reads_exp.stops(1:num_reads);
      paired_reads_exp.mates = paired_reads_exp.mates(:,1:num_reads/2);
    else
      paired_reads_exp.starts = []; paired_reads_exp.stops = []; paired_reads_exp.mates = [];
    end
  else
    % insert size set to median of empirical distribution
    ins_size = median(CFG.ins_sizes);
    max_iter = 1;
    cand1_starts = tidx(1:end-2*CFG.read_len-ins_size+1);
    cand1_stops = tidx(CFG.read_len-1+[1:end-2*CFG.read_len-ins_size+1]);
    cand2_starts = tidx(CFG.read_len+ins_size+[1:end-2*CFG.read_len-ins_size+1]);
    cand2_stops = tidx(2*CFG.read_len+ins_size-1+[1:end-2*CFG.read_len-ins_size+1]);
    paired_reads_exp.starts = [cand1_starts, cand2_starts];
    paired_reads_exp.stops = [cand1_stops, cand2_stops];
    paired_reads_exp.mates = [1:length(cand1_starts); length(cand1_starts)+1:length(cand1_starts)+length(cand2_starts)];
  end
  % connectivity matrix for expected pairs
  pair_mat_exp(:,:,t) = gen_paired_segments(segments, paired_reads_exp)./max_iter;
  pair_mat_exp(:,:,t) = pair_mat_exp(:,:,t)./pair_mat_len;
end
norm_pe = sum(sum(sum(pair_mat_exp)));
if norm_pe~=0 
  pair_mat_exp = pair_mat_exp./norm_pe;
end
assert(all(all(all(~isnan(pair_mat_exp)))));

%%% connectivity  matrix for observed pairs
pair_mat_obs = gen_paired_segments(segments, paired_reads_obs);
pair_mat_obs = pair_mat_obs./pair_mat_len;
