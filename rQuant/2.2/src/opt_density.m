function [profile_weights, obj] = opt_density(CFG, genes, profile_weights)
% OPT_DENSITY   Determines optimal parameters of density function.
%
%   [profile_weights, obj] = opt_profiles_smo(CFG, genes, profile_weights)
%
%   -- input --
%   CFG:             configuration struct
%   genes:           struct defining genes with start, stops, exons etc.
%   profile_weights: weights of profile functions (initialisation)
%
%   -- output --
%   profile_weights: weights of profile functions 
%   obj:             objective evaluated with optimal parameters
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


T = length([genes.transcripts]);       % number of transcripts
P = sum([genes.exonic_len]);           % number of positions
P_all = P;
F = CFG.num_plifs;                     % number of supporting points
N = size(CFG.transcript_len_ranges,1); % number of length bins

I = 0; % upper bound for number of introns
for g = 1:length(genes),
  introns = zeros(0,2);
  for t = 1:length(genes(g).transcripts),
    introns = [introns; genes(g).exons{t}(1:end-1,2)+1, genes(g).exons{t}(2:end,1)-1];
  end
  introns = unique(introns, 'rows');
  I = I + size(introns,1);
  clear introns;
end

profile_weights_adj = get_adj_bins(CFG);

%%%%% pre-processing %%%%%
coverage = zeros(P, 1);
exon_feat = sparse(P, T); % stores exon features from all transcripts in profile gene set
exon_feat_val = zeros(P, 1);
exon_feat_val_next = zeros(P, 1);
exon_feat_row = [1:P]';
exon_feat_col = zeros(P, 1);
intron_count = zeros(I, 1);
intron_mask = zeros(I, T);
mask = true(P, 1); 
ci = 0; ct = 0; cn = 0;
weights = zeros(1,T);
if CFG.VERBOSE>1, fprintf(1, 'Loading reads...\n'); tic; end
tmp_VERBOSE = CFG.VERBOSE;
CFG.VERBOSE = 0;
for g = 1:length(genes),
  if tmp_VERBOSE>1, fprintf(1, '%i\r', g); end
  try
    [tmp_coverage reads_ok tmp_introns] = get_coverage_per_read(CFG, genes(g), 1);
  catch
    reads_ok = 0;
  end
  assert(reads_ok==1);
  coverage(ci+[1:genes(g).exonic_len],1) = tmp_coverage;
  clear tmp_coverage;
  Tg = length(genes(g).transcripts);
  % initialisation of transcript weights to proportionate mean coverage
  weights(ct+[1:Tg]) = full(mean(coverage(ci+[1:genes(g).exonic_len]))/Tg*ones(1,Tg));
  for t = 1:length(genes(g).transcripts), % only works for tscp_len=1
    % exon features
    [exon_feat(ci+[1:genes(g).exonic_len], ct+1), tmp_exon_feat_val, tmp_exon_feat_val_next, feat_del_idx] = gen_exon_features(CFG, genes(g), t, 1);
    exon_feat_val(ci+[1:genes(g).exonic_len]) = tmp_exon_feat_val + F*ct; % transformation to linear indices
    exon_feat_val_next(ci+[1:genes(g).exonic_len]) = tmp_exon_feat_val_next + F*ct;
    exon_feat_col(ci+[1:genes(g).exonic_len]) = ct+1;
    genes(g).transcript_len_bin(t) = find(CFG.transcript_len_ranges(:,1) <= genes(g).transcript_length(t) & ...
                                          CFG.transcript_len_ranges(:,2) >= genes(g).transcript_length(t));
    clear tmp_exon_feat_val tmp_exon_feat_val_next;
  end
  % introns
  [tmp_intron_mask tmp_intron_count] = get_intron_data(genes(g), CFG, tmp_introns, g);
  intron_mask(cn+[1:length(tmp_intron_count)], ct+[1:Tg]) = tmp_intron_mask; 
  intron_count(cn+[1:length(tmp_intron_count)],1) = tmp_intron_count;
  % repeat mask
  fname = sprintf('%s%s_repeat', CFG.repeats_fn, genes(g).chr);
  if exist(sprintf('%s.pos', fname), 'file')
    [map.pos map.repeats] = interval_query(fname, {'repeats'}, [genes(g).start;genes(g).stop]);
    if ~isempty(map.pos)
      [tmp idx1 idx2] = intersect(map.pos, genes(g).eidx);
      assert(length(idx2)<=length(map.pos));
      mask(ci+idx2) = false;
    end
  end
  mask(ci+feat_del_idx) = false;
  ci = ci + genes(g).exonic_len;
  ct = ct + length(genes(g).transcripts);
  cn = cn + length(tmp_intron_count);
  clear tmp_introns tmp_intron_mask tmp_intron_count feat_del_idx;
end
assert(P==size(exon_feat,1));
% cut intron data to actual number of introns
if I < cn
  I = cn;
  assert(sum(sum(intron_mask(I+1:end,:)))==0);
  assert(sum(intron_count(I+1:end,1))==0);
  intron_mask = intron_mask(1:I,:);
  intron_count = intron_count(1:I,1);
end
CFG.VERBOSE = tmp_VERBOSE;
if CFG.VERBOSE>1, fprintf(1, 'Took %.1fs.\n', toc); end
tscp_len_bin = [genes.transcript_len_bin];
% find thetas that do not need to be optimised (located in the body of the profile function)
pw_nnz = get_included_thetas(CFG);

% subsample positions
if CFG.subsample,
  subsample_frac = min(CFG.subsample_frac, CFG.max_num_train_exm/sum(mask));
  tmp_P = round(P*subsample_frac);
  midx = find(mask);
  ridx = randperm(length(midx));
  mask(midx(ridx(tmp_P+1:end))) = false;
  clear midx ridx;
end

% exclude positions that are repetitive or subsampled or not representable by PLiFs
if any(~mask),
  subs_idx = find(mask);
  P = length(subs_idx);
  coverage = coverage(subs_idx, :);
  exon_feat = exon_feat(subs_idx, :);
  exon_feat_val = exon_feat_val(subs_idx, 1);
  exon_feat_val_next = exon_feat_val_next(subs_idx, 1);
  exon_feat_row = [1:P]';
  exon_feat_col = exon_feat_col(subs_idx, 1);
  if CFG.VERBOSE>0, fprintf('Subsampled from %i to %i positions\n', P_all, P); end
  clear P_old;
end

% adjust regularisation parameters
CFG.C_I = CFG.C_I*P/1000;
CFG.C_F = CFG.C_F*P/1000;
CFG.C_N = CFG.C_N*P/1000;


%%%%% optimisation %%%%%
eps = 1e-2;
C_w = [genes.transcript_length]';
%%% initialisation of variables
weights_old = zeros(1,T);
num_opt_steps = 2;
if nargin<3
  profile_weights = zeros(F, N);
  profile_weights(pw_nnz) = 1;
else
  profile_weights(~pw_nnz) = 0;
end
if CFG.VERBOSE>0 && CFG.load_profiles,
  fprintf(1, 'Using these profiles as initialisation\n')
  profile_weights
end
profile_weights_old = zeros(F, N);
norm_pw = ones(1, N);
fval = 1e100*ones(1, num_opt_steps); % objective values
fval_old = 0;
iter = 1;

if CFG.VERBOSE>0, fprintf('\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\t\tDelta norm\n'); end
while 1
  tic
  weights_old = weights;
  fval_old = fval;
  profile_weights_old = profile_weights;
  cnt = 1;
 
  %%%%% A. optimise transcript weights
  exon_mask = gen_exon_mask(profile_weights, tscp_len_bin, exon_feat, exon_feat_val, exon_feat_val_next, exon_feat_row, exon_feat_col);
  tmp_VERBOSE = CFG.VERBOSE;
  CFG.VERBOSE = 0;
  R_const = CFG.C_N*sum(sum((profile_weights(:,1:end-1)-profile_weights(:,2:end)).^2) + CFG.C_F*sum(sum((profile_weights(1:end-1,:)-profile_weights(2:end,:)).^2)));
  [weights, fval(cnt)] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, intron_mask, [], [], C_w, R_const, 1, weights, 'L1');
  %if ~(fval_old(end)-fval(cnt)>-1e-3)
  %  %fval_old(end)-fval(cnt)
  %end
  if cnt==2
    assert(fval(cnt-1)-fval(cnt)>-1e-3);
  end
  cnt = cnt + 1;
  CFG.VERBOSE = tmp_VERBOSE;
    
  %%%%% B. optimise profile weights
  profile_weights = reshape(profile_weights, 1, F*N);
  R_const = abs(weights)*C_w;
  if size(intron_mask,1)>0, R_const = R_const + CFG.C_I*sum((intron_mask*weights'-intron_count).^2); end
  tmp_VERBOSE = CFG.VERBOSE;
  CFG.VERBOSE = 0;
  [profile_weights, fval(cnt)] = opt_profiles_descent(CFG, profile_weights, tscp_len_bin, exon_feat, exon_feat_val, exon_feat_val_next, exon_feat_row, exon_feat_col, weights, coverage, [], R_const);
  % normalisation of weights
  profile_weights_unnorm = profile_weights;
  norm_pw = median(profile_weights(pw_nnz))*ones(1, N);
  profile_weights_unnorm = profile_weights;
  for n = 1:N,
    %norm_pw(n) = norm(profile_weights(:,n));
    profile_weights(:,n) = profile_weights(:,n)./norm_pw(n);
  end
  weights = weights.*norm_pw(tscp_len_bin);
  % recompute objective value
  exon_mask = gen_exon_mask(profile_weights, tscp_len_bin, exon_feat, exon_feat_val, exon_feat_val_next, exon_feat_row, exon_feat_col);
  p_adj_f = get_adj_bins(CFG);
  pw_nnz = get_included_thetas(CFG);
  pw_nnz1 = pw_nnz; pw_nnz1(F*(N-1)+1:end) = false;
  pw_nnz2 = pw_nnz; pw_nnz2(1:F) = false;
  fidx = find(pw_nnz1(1:F*(N-1)) & pw_nnz2(F+1:end));
  profile_weights = reshape(profile_weights, 1, F*N);
  R_const = abs(weights)*C_w;
  if size(intron_mask,1)>0, R_const = R_const + CFG.C_I*sum((intron_mask*weights'-intron_count).^2); end
  fval(cnt) = sum((exon_mask*weights'-coverage).^2) + R_const + CFG.C_N*sum((profile_weights(fidx)-profile_weights(fidx+F)).^2) + CFG.C_F*sum((profile_weights(pw_nnz)-profile_weights(p_adj_f(2,pw_nnz))).^2);
  profile_weights = reshape(profile_weights, F, N);
  if ~(fval(cnt-1)-fval(cnt)>-1e-3)
    %fval(cnt-1)-fval(cnt)
  end
  %assert(fval(cnt-1)-fval(cnt)>-1e-3);
  cnt = cnt + 1;
  CFG.VERBOSE = tmp_VERBOSE;
  
  %%%%% convergence criteria
  norm_weights = norm([weights_old, reshape(profile_weights_old,1,F*N)] - [weights, reshape(profile_weights,1,F*N)]);
  if fval_old(end)>=fval(end), sg = '-'; else sg = '+'; end
  if CFG.VERBOSE>1, fprintf(1, '%i\t%.3d\t%.5f\t%.3d\t%s\t%.1f\n', iter, fval(end), fval(end)/fval_old(end), norm_weights, sg, toc); end
  
  if norm(fval_old-fval)<eps || norm_weights<eps || iter >= CFG.max_iter
    break;
  end
  iter = iter + 1;
end
if CFG.VERBOSE>0, fprintf('Took %.1fs.\n', toc); end

obj = fval(end);
