function feats = append_pair_feats(chunks, feats, CFG)

% feats = append_pair_feats(chunks, feats, CFG)
%
% 
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009-2011

assert(issorted(chunks, 'rows'));

tic
L = sum(chunks(:,3)-chunks(:,2)+1);
% features correspond to a horizontal concatenation (of length L) of feature
% blocks of all individual chunks

% number of read pairs aligned left and right of the given position (with
% alignment parameters indicating that these are paired reads)
pair_cover_feat = zeros(2, L);

strands = '+-';
offset = 0;

tic
for c=1:size(chunks,1),
  region.start   = chunks(c,2);
  region.stop    = chunks(c,3);
  region.chr_num = chunks(c,1);
  region.chr = CFG.chr_names{chunks(c,1)};
  for s=1:length(strands),
    region.strand = strands(s);
    pair_cover = get_pair_features_from_bam(CFG, region);

    chunk_idx = offset+1 : offset+region.stop-region.start+1;
    assert(all(pair_cover_feat(s,chunk_idx)  == 0));
    pair_cover_feat(s,chunk_idx) = double(pair_cover);
  end
  % advance offset to point to the start of the next feature chunk
  offset = offset + region.stop - region.start + 1;
%  fprintf('  processed %i (of %i) chunks, %2.1f%% in %.1f sec\r', c, ...
%          size(chunks,1), 100*c/size(chunks,1), toc)
end
%fprintf('  processed %i (of %i) chunks, %2.1f%% in %.1f sec\n', c, ...
%        size(chunks,1), 100*c/size(chunks,1), toc)

% first feature is total pair span, second feature is a summation of
% all coverage/span features indicating gene presence, third feature
% indicates the length of completely uncovered blocks 
pair_feat = zeros(3, size(pair_cover_feat,2));
pair_feat(1,:) = max(pair_cover_feat);
% add up exon coverage, intron span and mate-pair span features
pair_feat(2,:) = pair_feat(1,:) + feats(1,:) + feats(4,:);
low_cover_blocks = find_blocks(pair_feat(2,:) < CFG.low_cover_threshold);
l = low_cover_blocks(2,:) - low_cover_blocks(1,:) + 1;
for b=1:size(low_cover_blocks,2),
  idx = low_cover_blocks(1,b):low_cover_blocks(2,b);
  pair_feat(3, idx) = log2(l(b) / mean(pair_feat(2,idx)+1));
end

feats = [feats; pair_feat];
assert(size(feats,2) == L);
% assert that the total coverage feature is indeed at least as high as the
% sum of all relevant coverage / span features
assert(all(feats(end-1,:) >= feats(1,:)+feats(4,:)+feats(end-2,:)));

fprintf('  appended read-pair features to %i chunks in %.1f sec\n', size(chunks,1), toc);

% eof