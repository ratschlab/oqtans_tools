function [profile_weights, obj, fval] = opt_profiles_descent(CFG, profile_weights, tscp_len_bin, exon_feat, exon_feat_val, exon_feat_val_next, exon_feat_row, exon_feat_col, weights, coverage, seq_coeff, R_const)
% OPT_PROFILES_DESCENT   Determines the optimal profile functions.
%
%   [profile_weights, obj, fval] = opt_profiles_descent(CFG, profile_weights, tscp_len_bin, exon_feat, exon_feat_val, exon_feat_val_next, exon_feat_row, exon_feat_col, weights, coverage, seq_coeff, R_const)
%
%   -- input --
%   CFG:                configuration struct
%   profile_weights:    weights of profile functions
%   tscp_len_bin:       vector of length bins for each transcript
%   exon_feat:          P x num_bins matrix of features for P exonic positions
%   exon_feat_val:      vector of indices to supporting points f
%   exon_feat_val_next: vector of indices to supporting points f+1
%   exon_feat_row:      vector of row indices for sparse P x T matrix (position index)
%   exon_feat_col:      vector column indices for sparse P x T matrix (transcript index)
%   weights:            weights of transcripts
%   coverage:           vector of observed exon coverage
%   seq_coeff:          vector of sequence correction
%   R_const:            constant residue
%
%   -- output --
%   profile_weights:    weights of profile functions
%   obj:                objective value at optimum
%   fval:               objective value at each step
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


F = CFG.num_plifs;                     % number of supporting points
N = size(CFG.transcript_len_ranges,1); % number of length bins
T = length(weights);
tp_idx = zeros(1,F*T);
for n = 1:F,
  tp_idx([1:F:length(tp_idx)]+n-1) = n+[0:(F*T+F):F*T*T];
end

max_iter = 1;

profile_weights = reshape(profile_weights, 1, F*N);
profile_weights_old = zeros(1, F*N);

pw_nnz = get_included_thetas(CFG);
pw_nnz = reshape(pw_nnz, 1, F*N);
pidx = find(pw_nnz);

% adjacent supporting points
p_adj_f = get_adj_bins(CFG);

% adjacent transcript length bins
%p_adj_n = zeros(2, F*N);
%p_adj_n(2,:) = (F+1):F*(N+1);
%p_adj_n(2,F*N-F+1:F*N) = p_adj_n(2,F*(N-2)+1:F*(N-1)); 

fval = 1e100*ones(1, length(pidx));
fval_old = zeros(1, length(pidx));

if CFG.VERBOSE>0, fprintf(1, '\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tNorm diff\n'); end
iter = 1;
while 1
  cnt = 1;
  profile_weights_old = profile_weights;
  fval_old = fval;
  for p = pidx,
    %fprintf('%i\r', p);
    f1 = mod(p,F); if f1==0, f1 = F; end
    n1 = ceil(p/F);
    % entries corresponding to theta_f_n
    if f1==F
      idx_w_th1 = find(mod(exon_feat_val,F)==0 & tscp_len_bin(exon_feat_col)'==n1);
    else
      idx_w_th1 = find(mod(exon_feat_val,F)==f1 & tscp_len_bin(exon_feat_col)'==n1);
    end
    % entries corresponding to theta_f-1_n
    if f1==1
      idx_w_th2 = zeros(0,1);
    else
      idx_w_th2 = find(mod(exon_feat_val,F)==f1-1 & tscp_len_bin(exon_feat_col)'==n1);
    end
    % entries corresponding to theta_not(f/f-1,n)
    idx_wo_th = setdiff(1:length(exon_feat_val), [idx_w_th1; idx_w_th2])';
    exon_feat_ones = sparse(exon_feat_row, exon_feat_col, 1, length(exon_feat_col), T);
    [exon_mask, exon_mask_part1, exon_mask_part2] = gen_exon_mask(reshape(profile_weights, F, N), tscp_len_bin, exon_feat, exon_feat_val, exon_feat_val_next, exon_feat_row, exon_feat_col);
    %%% Rth1: residue for theta_f_n
    Rth1 = (exon_feat_ones(idx_w_th1,:)-exon_feat(idx_w_th1,:))*weights';
    %%% Rth2: residue for theta_f-1_n
    Rth2 = exon_feat(idx_w_th2,:)*weights';
    %%% R1: residue for theta_f_n independent variables
    R1 = exon_mask_part2(idx_w_th1,:)*weights' - coverage(idx_w_th1);
    %%% R2: residue for theta_f-1_n independent variables
    R2 = exon_mask_part1(idx_w_th2,:)*weights' - coverage(idx_w_th2);
    %%% R3: residue for theta_f/f-1_n independent variables
    R3 = exon_mask(idx_wo_th,:)*weights' - coverage(idx_wo_th);
    clear exon_mask exon_mask_part1 exon_mask_part2;
    %%% R4: residue for coupling transcript length bins
    R4 = 0;
    if n1<N && pw_nnz(f1+n1*F), R4 = R4 + profile_weights(f1+n1*F)^2; end
    if n1>1 && pw_nnz(f1+(n1-2)*F), R4 = R4 + profile_weights(f1+(n1-2)*F)^2; end
    for f = 1:F,
      for n = 1:N-1,
        if (f~=f1 || (f==f1 && n~=n1 && n~=n1-1)) && pw_nnz(f+(n-1)*F) && pw_nnz(f+n*F)
          R4 = R4 + (profile_weights(f+(n-1)*F)-profile_weights(f+n*F))^2;
        end
      end
    end
    %%% R5: residue for coupling supporting points
    R5 = 0;
    if f1<F, R5 = R5 + profile_weights(p_adj_f(2,f1+(n1-1)*F))^2; end
    if f1>1, R5 = R5 + profile_weights(p_adj_f(1,f1+(n1-1)*F))^2; end
    for n = 1:N,
      for f = 1:F-1,
        if (n~=n1 || (n==n1 && f~=f1 && p_adj_f(2,f+(n-1)*F)~=(f1+(n-1)*F))) && pw_nnz(f+(n-1)*F)
          %ii =[ii, [f+(n-1)*F; p_adj_f(2,f+(n-1)*F)]];
          R5 = R5 + (profile_weights(f+(n-1)*F)-profile_weights(p_adj_f(2,f+(n-1)*F)))^2;
        end
      end
    end
    %%% S1: residue of quadratic term
    S1 = sum(Rth1.^2) + sum(Rth2.^2);
    % coupling constraints
    if n1<N && pw_nnz(f1+n1*F), S1 = S1 + CFG.C_N; end
    if n1>1 && pw_nnz(f1+(n1-2)*F), S1 = S1 + CFG.C_N; end
    if f1<F, S1 = S1 + CFG.C_F; end
    if f1>1, S1 = S1 + CFG.C_F; end
    assert(S1>0); % condition for minimum (2nd derivative > 0)
    %%% S2: residue of linear term
    S2 = 2*sum(Rth1'*R1) + 2*sum(Rth2'*R2);
    % coupling constraints
    if n1<N && pw_nnz(f1+n1*F), S2 = S2 - 2*CFG.C_N*profile_weights(f1+n1*F); end
    if n1>1 && pw_nnz(f1+(n1-2)*F), S2 = S2 - 2*CFG.C_N*profile_weights(f1+(n1-2)*F); end
    if f1<F, S2 = S2 - 2*CFG.C_F*profile_weights(p_adj_f(2,f1+(n1-1)*F)); end
    if f1>1, S2 = S2 - 2*CFG.C_F*profile_weights(p_adj_f(1,f1+(n1-1)*F)); end
    %%% S3: constant term
    S3 = sum(R1.^2) + sum(R2.^2) + sum(R3.^2) + CFG.C_N*R4 + CFG.C_F*R5 + R_const;
    %%% calculation and clipping of theta
    th_new = -0.5*S2/S1;
    if th_new < 0
      profile_weights(p) = 0.0;
    else
      profile_weights(p) = th_new;
    end
    fval(cnt) = quad_fun(profile_weights(p), S1, S2, S3);
    tmp_pw = reshape(profile_weights, F, N);
    exon_mask = gen_exon_mask(reshape(profile_weights, F, N), tscp_len_bin, exon_feat, exon_feat_val, exon_feat_val_next, exon_feat_row, exon_feat_col);
    pw_nnz1 = pw_nnz; pw_nnz1(F*(N-1)+1:end) = false;
    pw_nnz2 = pw_nnz; pw_nnz2(1:F) = false;
    fidx = find(pw_nnz1(1:F*(N-1)) & pw_nnz2(F+1:end));
    obj_alt = sum((exon_mask*weights'-coverage).^2) + R_const + CFG.C_N*sum((profile_weights(fidx)-profile_weights(fidx+F)).^2) + CFG.C_F*sum((profile_weights(pw_nnz)-profile_weights(p_adj_f(2,pw_nnz))).^2);
    if ~(abs(fval(cnt)-obj_alt)<1e-3) % objective should be indentical to not-expanded objective
      if CFG.VERBOSE>1, fprintf(1, 'objectives differ %.6f (theta %i)\n', abs(fval(cnt)-obj_alt), p); end
    end
    cnt = cnt + 1;
  end
  if CFG.VERBOSE>1, fprintf(1, '%i\t%.5d\t%.5d\n', iter, fval(end), norm(profile_weights_old-profile_weights)); end
  if norm(fval_old-fval)<1e-5 || norm(profile_weights_old-profile_weights)<1e-5 || iter>=max_iter,
    break;
  end
  iter = iter + 1;
end
if 0%mean(fval(1:end-1)-fval(2:end)>-1e-2)<0.95
  fprintf(1, 'objective non-decreasing\n');
end
assert(all(fval(1:end-1)-fval(2:end)>-1e-3)); % objective should decrease at every step
if CFG.VERBOSE>0, fprintf(1, 'Took %.1fs.\n', toc); end

profile_weights = reshape(profile_weights, F, N);
obj = fval(end);