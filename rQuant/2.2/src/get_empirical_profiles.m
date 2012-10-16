function profile_weights = get_empirical_profiles(CFG, genes)
% GET_EMPIRICAL_PROFILES   Estimates profiles empirically from read data.
%
%   profile_weights = get_empirical_profiles(CFG, genes)
%
%   -- input --
%   CFG:             configuration struct
%   gene:            struct defining a gene with start, stops, exons etc.
%
%   -- output --
%   profile_weights: empirically estimated profiles
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


N = size(CFG.transcript_len_ranges,1);
profile_weights = zeros(CFG.num_plifs, N);
num_exm = zeros(1, N); t = 1;
pw_nnz = get_included_thetas(CFG);
lmt = get_limits(CFG.max_side_len, CFG.num_plifs/2);
for g = 1:length(genes),
  if length(genes(g).transcripts)>1, continue; end
  [tmp_coverage reads_ok] = get_coverage_per_read(CFG, genes(g), 1);
  assert(reads_ok==1);
  if genes(g).exonic_len>=CFG.num_plifs
    fidx = find(lmt<=genes(g).exonic_len/2, 1, 'last');
    take_idx = lmt(1:fidx);
    profile_weights(1:length(take_idx),genes(g).transcript_len_bin(t)) = profile_weights(1:length(take_idx),genes(g).transcript_len_bin(t)) + full(tmp_coverage(take_idx));
    take_idx = genes(g).exonic_len-lmt(fidx:-1:1)+1;
    profile_weights(end-length(take_idx)+1:end, genes(g).transcript_len_bin(t)) = profile_weights(end-length(take_idx)+1:end, genes(g).transcript_len_bin(t)) + full(tmp_coverage(take_idx));
    num_exm(genes(g).transcript_len_bin(t)) = num_exm(genes(g).transcript_len_bin(t)) + 1;
  end
end
for n = 1:N,
  profile_weights(:,n) = profile_weights(:,n)./num_exm(n);
end
norm_pw = median(profile_weights(pw_nnz))*ones(1, N);
%norm_pw = max(profile_weights,[],1);
for n = 1:N,
  profile_weights(:,n) = profile_weights(:,n)./norm_pw(n);
end