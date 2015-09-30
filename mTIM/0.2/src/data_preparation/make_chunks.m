function chunks = make_chunks(CFG)

% chunks = make_chunks(CFG)
%
% Subdivides chromosomes into chunks splitting only in intergenic regions
% that are not spanned by any read. The returned chunks have four
% columns: chromosome number, start position, end position, and chunk id.
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009-2011

USE_READS = 0 

MIN_IGE_BLOCK_SIZE = 100;
MAX_IGE_BLOCK_SIZE = 1e6; % not longer than 1mio nt

MIN_EXON_COVER = 5;
MIN_INTRON_COVER = 5;

genic_map = cell(1,CFG.num_chr);
for c=1:CFG.num_chr,
  genic_map{c} = zeros(1,CFG.chr_lens(c), 'uint8');
  exon_map{c}  = zeros(1,CFG.chr_lens(c), 'uint8');
end

% mask positions as genic that are labeled as such
if ~exist(CFG.label_fn, 'file')
   load([CFG.label_fn '.mat'], '-mat', 'label_map', 'LABELS');
else
   load(CFG.label_fn, 'label_map', 'LABELS');
end
for c=1:CFG.num_chr,
  idx = find(label_map{c} ~= LABELS.intergenic);
  genic_map{c}(idx) = 1;
  idx = find(label_map{c} == LABELS.exon_W | label_map{c} == LABELS.exon_C);
  exon_map{c}(idx) = 1;
end
clear LABELS label_map


if USE_READS,
  % mask positions as genic if they are covered by exon or intron reads
  % generate global chromosomal maps here
  strands = '+-';
  exon_read_map   = cell(1, CFG.num_chr);
  intron_read_map = cell(1, CFG.num_chr);

  tic
  for c=1:CFG.num_chr,
    exon_map{c}   = zeros(1, CFG.chr_lens(c), 'uint8');
    intron_map{c} = zeros(1, CFG.chr_lens(c), 'uint8');

    region.start   = 1;
    region.stop    = CFG.chr_lens(c);
    region.chr_num = c;
    for s=1:length(strands),
      region.strand = strands(s);
      [exon_cover, intron_span] = get_read_features_from_bam(CFG, region);
      % here we simply add up thr strand-specific read counts
      exon_map{c}   = exon_map{c}   + uint8(exon_cover);
      intron_map{c} = intron_map{c} + uint8(intron_span);
    end
    fprintf('  generated read coverage map for chromosome %s (%.1f sec)\r', ...
            CFG.chr_names{c}, toc)
  end 
  fprintf('  generated read coverage map for chromosome %s (%.1f sec)\n', ...
          CFG.chr_names{c}, toc)

  for c=1:CFG.num_chr,
    idx = find(exon_map{c} >= MIN_EXON_COVER);
    genic_map{c}(idx) = 1;
    idx = find(intron_map{c} >= MIN_INTRON_COVER);
    genic_map{c}(idx) = 1;
  end
end


% identify intergenic blocks
chunks = zeros(0, 4);
prev_id = 0;
for c=1:CFG.num_chr,
  fprintf('  making chunks for %s', CFG.chr_names{c});
  ige_blocks = find_blocks(genic_map{c} == 0);
  % do not split too small blocks
  ige_blocks(:,ige_blocks(2,:)-ige_blocks(1,:) < MIN_IGE_BLOCK_SIZE) = [];
  
  %split_cands = [0, round((ige_blocks(1,:) + ige_blocks(2,:)) / 2), CFG.chr_lens(c)];
  %splits = [1, round((ige_blocks(1,:) + ige_blocks(2,:)) / 2), CFG.chr_lens(c)];
  
  %{
  ATTENTION: FUNCTIONALITY WAS REMOVED
  - this piece of code checks if splits contain at least some
    exonic part 
  - however, some too large blocks turned out to be intergenic and
    have to be splitted as well

  is_good = zeros(size(splits));
  for s=2:length(splits)-1,
    if any(exon_map{c}(splits(s-1):splits(s))) ...
          && any(exon_map{c}(splits(s):splits(s+1))),
      is_good(s) = 1;
    end
  end
  splits = splits(is_good==1);
%}  
  
  split_cands = [0, round((ige_blocks(1,:) + ige_blocks(2,:)) / 2), CFG.chr_lens(c)];
  
  % split also too long intergenic blocks
  % (e.g. in drosophila such a block could be >2.9mio nt)
  splits = [split_cands(1)];
  for i=1:length(split_cands)-1,
    len = split_cands(i+1)-split_cands(i);
    if (len<=MAX_IGE_BLOCK_SIZE),
      % size is okay
      splits = [splits, split_cands(i+1)];
    else
      fprintf('\nSplit-length is too big: %i \n',len);
      % the size of the current chunk is too big
      % therefore split it into n times MAX_IGE_BLOCK_SIZE chunks
      val = split_cands(i) + MAX_IGE_BLOCK_SIZE ;
      mul = floor(len/MAX_IGE_BLOCK_SIZE);
      fprintf('  val=%i   mul=%i   #splits=%i\n',val,mul,length(splits));
      for j=1:mul,
        splits = [splits, val];
        val = val + MAX_IGE_BLOCK_SIZE;
      end
      splits = [splits, split_cands(i+1)];
      fprintf('  #splits=%i\n',length(splits));
      % now, the other way around: if the last split
      % is too small merge it with the previous chunk
      if ((splits(end)-splits(end-1))<(10*MIN_IGE_BLOCK_SIZE)),
        splits = [splits(1:end-2), splits(end)];
        fprintf('  merging last split   #splits=%i\n',length(splits));
      end
      fprintf(' %i, ',splits);
      fprintf(' \nfinished.\n');
    end
  end

  chr = c*ones(length(splits)-1,1);
  b = [splits(1:end-1)'+1];
  e = [splits(2:end)'];
  assert(all(b < e));
  assert(all(e(1:end-1)+1 == b(2:end)));
  ids = [prev_id+1:prev_id+length(splits)-1]';
  prev_id = prev_id+length(splits)-1;
  chunks = [chunks; ...
           [chr, b, e, ids]];
  fprintf(' (%i)...\n', length(splits));
end
assert(length(unique(chunks(:,4))) == size(chunks,1));


% find chunks that are still >MAX_IGE_BLOCK_SIZE and report
lens = abs(chunks(:,2)-chunks(:,3))+1;
fprintf('Max/min chunk length:  %i/%i.\n', max(lens),min(lens));

% eof
