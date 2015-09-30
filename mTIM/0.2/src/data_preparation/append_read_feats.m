function feats = append_read_feats(chunks, feats, CFG)

% feats = append_read_feats(chunks, feats, CFG)
%
% 
%
% written by Georg Zeller, Pramod Mudrakarta & Andre Kahles, MPI Tuebingen, Germany, 2009-2011

%assert(issorted(chunks, 'rows'));
tic

L = sum(chunks(:,3)-chunks(:,2)+1);
% features correspond to a horizontal concatenation (of length L) of feature
% blocks of all individual chunks

% number of reads aligned at a given position
exon_cover_feat   = zeros(2, L);
% number of "introns" from spliced reads spanning a given position
intron_span_feat  = zeros(2, L);
% score derived from spliced-read alignments indicating a splice junction
% at a given position (score is proportional to the min. number of read
% positions aligned to either end)
spliced_read_feat = zeros(4, L);

strands = '+-';
offset = 0;

for c=1:size(chunks,1),
  region.start   = chunks(c,2);
  region.stop    = chunks(c,3);
  region.chr_num = chunks(c,1);
  region.chr = CFG.chr_names{chunks(c,1)};
  for s=1:length(strands),
    region.strand = strands(s);
    [exon_cover intron_span intron_list] = get_read_features_from_bam(CFG, region);

    chunk_idx = offset+1 : offset+region.stop-region.start+1;
    assert(all(exon_cover_feat(s,chunk_idx)  == 0));
    assert(all(intron_span_feat(s,chunk_idx) == 0));
    exon_cover_feat(s,chunk_idx)  = double(exon_cover);
    intron_span_feat(s,chunk_idx) = double(intron_span);

    assert(all(region.start <= intron_list(:,1) & intron_list(:,2) <= region.stop));
    idx_b = intron_list(:,1)' - region.start + 1 + offset;
    idx_e = intron_list(:,2)' - region.start + 1 + offset;
    assert(all(idx_e - idx_b > 0))
    % instead of this score (derived from the length of the length of the
    % spliced alignments in either neighboring exon) one could also use
    % the count of supporting introns (intron_list(:,4)') 
    score = intron_list(:,5)';
    spliced_read_feat(s*2-1,idx_b) = double(score);
    spliced_read_feat(s*2,  idx_e) = double(score);
  end
  % advance offset to point to the start of the next feature chunk
  offset = offset + region.stop - region.start + 1;
%  fprintf('  processed %i (of %i) chunks, %2.1f%% in %.1f sec\r', c, ...
%          size(chunks,1), 100*c/size(chunks,1), toc)
end
%fprintf('  processed %i (of %i) chunks, %2.1f%% in %.1f sec\n', c, ...
%        size(chunks,1), 100*c/size(chunks,1), toc)

% transform exon coverage features into one strand-independent coverage and
% a strand-discriminating difference feature
exon_feat = [max(exon_cover_feat); exon_cover_feat(1,:)-exon_cover_feat(2,:)];
% additionally, compute a coverage gradient feature indicating the magnitude
% of local change between lambda positions up- and downstream respectively
lambda = 2;
% computing this gradient w.r.t. to log2-coverage puts more weight on
% changes from/to 0-coerage
cover_gradient_feat = log2(1+exon_feat(1,2*lambda+1:end)) - log2(1+exon_feat(1,1:end-2*lambda));
cover_gradient_feat = [zeros(1,lambda), cover_gradient_feat, zeros(1,lambda)];
assert(size(cover_gradient_feat,2) == size(exon_feat,2));
exon_feat = [exon_feat; cover_gradient_feat];
% transform intron span features into one strand-independent span and a
% strand-discriminating difference feature
intron_feat = [max(intron_span_feat); intron_span_feat(1,:)-intron_span_feat(2,:)];

feats = [feats; exon_feat; intron_feat; spliced_read_feat];
assert(size(feats,2) == L);

fprintf('  appended read features to %i chunks in %.1f sec\n', size(chunks,1), toc);

% eof
