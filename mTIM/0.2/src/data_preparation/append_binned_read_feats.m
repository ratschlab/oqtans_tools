function feats = append_binned_read_feats(chunks, feats, CFG)
% feats = append_binned_read_feats(chunks, feats, CFG)
%
% Most basic version:
% - only handles intron spans (binary)
% - 3 bins, manually defined ranges
%
% written by Nico Görnitz, 03-2012
% based on append_read_feats 
% by Georg Zeller, Pramod Mudrakarta & Andre Kahles, MPI Tuebingen, Germany, 2009-2011

%assert(issorted(chunks, 'rows'));
tic

L = sum(chunks(:,3)-chunks(:,2)+1);
% features correspond to a horizontal concatenation (of length L) of feature
% blocks of all individual chunks

strands = '+-';
offset = 0;

bak = CFG.read_intron_span_min_score;

% for future extensions sort the levels
% array to be in ascending order
levels = [0,10,25];
levels = [0,20,35]; % new levels for d.melanogaster 130107
levels = sort(levels,'ascend');

intron_feat = [];
bins = {};

% get intron spans for different levels of 
% span_min_score (assume that levels is 
% in ascending order).
for i=1:length(levels),
    CFG.read_intron_span_min_score = levels(i);
    bins{i} = getIntronFeat(CFG, strands, offset, chunks, L);
end

% Subtract the higher levels from the lower to
% get a clean representation.
for i=length(levels)-1:-1:1,
    for j=1:i,
        bins{j} = bins{j}-bins{i+1};
    end
end

% Add the features in ascending order.
for i=1:length(levels),
    % transform intron span features into one strand-independent span and a
    % strand-discriminating difference feature
    intron_feat = [intron_feat; max(bins{i})];
end

feats = [feats; intron_feat];
assert(size(feats,2) == L);

CFG.read_intron_span_min_score = bak;
fprintf('  appended binned read features to %i chunks in %.1f sec\n', size(chunks,1), toc);




function [intron_span_feat] = getIntronFeat(CFG, strands, offset, chunks, L)
% return the number of "introns" from 
% spliced reads spanning a given position
intron_span_feat  = zeros(2, L);
for c=1:size(chunks,1),
  region.start   = chunks(c,2);
  region.stop    = chunks(c,3);
  region.chr_num = chunks(c,1);
  region.chr = CFG.chr_names{chunks(c,1)};
  for s=1:length(strands),
    region.strand = strands(s);
    [exon_cover intron_span intron_list] = get_read_features_from_bam(CFG, region);

    chunk_idx = offset+1 : offset+region.stop-region.start+1;
    assert(all(intron_span_feat(s,chunk_idx) == 0));
    intron_span_feat(s,chunk_idx) = double(intron_span);
  end
  % advance offset to point to the start of the next feature chunk
  offset = offset + region.stop - region.start + 1;
end


% eof
