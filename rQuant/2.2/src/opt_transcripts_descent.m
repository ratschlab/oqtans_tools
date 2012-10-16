function [weights, obj, fval] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, intron_mask, paired_exp, paired_obs, C_w, R_const, max_iter, weights0, reg)
% OPT_TRANSCRIPTS_DESCENT   Determines the optimal transcript weights.
%
%   [weights, obj, fval] = opt_transcripts_descent(CFG, coverage, exon_mask, intron_count, intron_mask, paired_exp, paired_obs, C_w, R_const, max_iter, weights0, reg)
%
%   -- input --
%   CFG:          configuration struct
%   coverage:     vector of observed exon coverage 
%   exon_mask:    PxT matrix modelling 'importance' of a position within a transcript
%   intron_count: vector of observed intron confirmation
%   intron_mask:  IxT matrix defining whether an intron belongs to a particular transcript
%   paired_exp:   matrix (#segments x #segments x T) of connectivity by expected paired-end reads
%   paired_obs:   matrix (#segments x #segments) of connectivity by observed paired-end reads
%   C_w:          regularisation parameter per transcript (T x 1)
%   R_const:      constant residue
%   max_iter:     maximal number of iterations (optional)
%   weights0:     initialisation values of the weights (optional)
%   reg:          regularisation method (optional)
%
%   -- output --
%   weights:      weights of transcripts
%   obj:          objective value at optimum
%   fval:         objective value at each step
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


T = size(exon_mask,2); % number of transcripts
I = size(intron_mask,1); % number of introns
PE = size(paired_exp,1); % number of segments

exon_count = sum(coverage,2);

if nargin<9
  R_const = 0;
end

if nargin<10
  max_iter = 1e100;
end

if nargin<11
  weights = full(mean(coverage)/T*ones(1,T));
else
  weights = weights0;
end
weights_old = zeros(1,T);
fval = 1e100*ones(1,T);
fval_old = zeros(1,T);

if nargin<12
  reg = 'L1';
end

LB = 0.0;
UB = full(mean(coverage));

cnt = 0;
if CFG.VERBOSE>0, fprintf(1, '\nStarting optimising...\n'); tic; end
if CFG.VERBOSE>1, fprintf(1, 'Itn\tObjective\tNorm diff\n'); end
if T==1
  t = 1;
  % exon loss
  RE = -exon_count;
  S1 = sum(exon_mask(:,t).^2);
  S2 = 2*sum(exon_mask(:,t)'*RE);
  S3 = sum(RE.^2) + R_const;
  % intron loss
  if I>0
    RI = -intron_count;
    S1 = S1 + CFG.C_I*sum(intron_mask(:,t).^2);
    S2 = S2 + 2*CFG.C_I*sum(intron_mask(:,t)'*RI);
    S3 = S3 + CFG.C_I*sum(RI.^2);
  end
  % paired loss
  if PE>0
    RPE = -paired_obs;
    S1 = S1 + CFG.C_PE*sum(sum((squeeze(paired_exp(:,:,t))).^2));
    S2 = S2 + 2*CFG.C_PE*sum(sum(squeeze(paired_exp(:,:,t)).*RPE));
    S3 = S3 + CFG.C_PE*sum(sum(RPE.^2));
  end
  % regularisation
  switch reg
   case 'L1'
    S2 = S2 + C_w(t);
   case 'L2'
    S1 = S1 + C_w(t);
  end
  w_new = -0.5*S2/S1;
  % clipping of w_t
  if w_new < 0
    weights(t) = LB;
  else
    weights(t) = w_new;
  end
  fval(t) = quad_fun(weights(t), S1, S2, S3);
  obj_alt = sum((exon_mask*weights'-coverage).^2) + R_const;
  if I>0, obj_alt = obj_alt + CFG.C_I*sum((intron_mask*weights'-intron_count).^2); end
  if PE>0
    paired_exp_wsum = squeeze(paired_exp(:,:,t))*weights(t);
    obj_alt = obj_alt + CFG.C_PE*sum(sum((paired_exp_wsum-paired_obs).^2));
  end
  switch reg
   case 'L1'
    obj_alt = obj_alt + abs(weights*C_w);
   case 'L2'
    obj_alt = obj_alt + weights.^2*C_w;
  end
  if ~(abs(fval(t)-obj_alt)<1e-3) % objective should be indentical to not-expanded objective
    cnt = cnt + 1;
    if CFG.VERBOSE>1, fprintf(1, 'objectives differ %.6f (tscp %i)\n', abs(fval(t)-obj_alt), t); end
  end
else
  iter = 1;
  while 1
    weights_old = weights;
    fval_old = fval;
    for t = 1:T,
      idx_wo_t = setdiff(1:T,t);
      % exon loss
      RE = exon_mask(:,idx_wo_t)*weights(idx_wo_t)'-exon_count;
      S1 = sum(exon_mask(:,t).^2);
      S2 = 2*sum(exon_mask(:,t)'*RE);
      S3 = sum(RE.^2) + R_const;
      % intron loss
      if I>0
        RI = intron_mask(:,idx_wo_t)*weights(idx_wo_t)'-intron_count;
        S1 = S1 + CFG.C_I*sum(intron_mask(:,t).^2);
        S2 = S2 + 2*CFG.C_I*sum(intron_mask(:,t)'*RI);
        S3 = S3 + CFG.C_I*sum(RI.^2);
      end
      % paired loss
      if PE>0
        RPE = -paired_obs;
        for tt = idx_wo_t,
          RPE = RPE + squeeze(paired_exp(:,:,tt))*weights(tt);
        end
        S1 = S1 + CFG.C_PE*sum(sum((squeeze(paired_exp(:,:,t))).^2));
        S2 = S2 + 2*CFG.C_PE*sum(sum(squeeze(paired_exp(:,:,t)).*RPE));
        S3 = S3 + CFG.C_PE*sum(sum(RPE.^2));
      end
      % regularisation
      switch reg
       case 'L1'
        S2 = S2 + C_w(t);
        S3 = S3 + abs(weights(idx_wo_t))*C_w(idx_wo_t);
       case 'L2'
        S1 = S1 + C_w(t);
        S3 = S3 + weights(idx_wo_t).^2*C_w(idx_wo_t);
      end
      w_new = -0.5*S2/S1;
      % clipping of w_t
      if w_new < 0
        weights(t) = LB;
      else
        weights(t) = w_new;
      end
      fval(t) = quad_fun(weights(t), S1, S2, S3);
      obj_alt = sum((exon_mask*weights'-coverage).^2) + R_const;
      if I>0, obj_alt = obj_alt + CFG.C_I*sum((intron_mask*weights'-intron_count).^2); end
      if PE>0,
        paired_exp_wsum = zeros(PE, PE);
        for tt = 1:T,
          paired_exp_wsum = paired_exp_wsum + squeeze(paired_exp(:,:,tt))*weights(tt);
        end
        obj_alt = obj_alt + CFG.C_PE*sum(sum((paired_exp_wsum-paired_obs).^2));
      end
      switch reg
       case 'L1'
        obj_alt = obj_alt + abs(weights*C_w);
       case 'L2'
        obj_alt = obj_alt + weights.^2*C_w;
      end
      if ~(abs(fval(t)-obj_alt)<1e-3) % objective should be indentical to not-expanded objective
        cnt = cnt + 1;
        if CFG.VERBOSE>1, fprintf(1, 'objectives differ %.6f (tscp %i)\n', abs(fval(t)-obj_alt), t); end
      end
    end
    if CFG.VERBOSE>1, fprintf(1, '%i\t%.5d\t%.5d\n', iter, fval(end), norm(weights_old-weights)); end
    if norm(fval_old-fval)<1e-5 || norm(weights_old-weights)<1e-5 || iter>=max_iter,
      break;
    end
    iter = iter + 1;
  end
end
if CFG.VERBOSE>0 && cnt>0, fprintf(1, 'objectives differ for %i transcripts\n', cnt); end
if mean(fval(1:end-1)-fval(2:end)>-1e-2)<0.95
  fprintf(1, 'objective non-decreasing\n');
end
%assert(all(fval(1:end-1)-fval(2:end)>-1e-2));
if CFG.VERBOSE>0, fprintf(1, 'Took %.1fs.\n', toc); end

obj = fval(end);