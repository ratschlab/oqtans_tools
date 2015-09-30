function chunks = make_genome_chunks(CFG)

chunk_size = 5000000; % todo
overlap    =  500000; % todo

chunks = zeros(0, 4);
prev_id = 0;
for c=1:CFG.num_chr,
  fprintf('  making genome chunks for chromosome %s', CFG.chr_names{c});
  N = ceil(CFG.chr_lens(c) / chunk_size);
  splits = round(linspace(1, CFG.chr_lens(c), N+1));
  splits = splits(2:end-1);
  
  chr = c*ones(length(splits)+1,1);
  b = [1; splits'+1-overlap/2];
  e = [splits'+overlap/2; CFG.chr_lens(c)];
  assert(all(b < e));
  assert(all(e(1:end-1) > b(2:end)));
  ids = [prev_id+1:prev_id+length(splits)+1]';
  prev_id = prev_id+length(splits)+1;
  chunks = [chunks; ...
           [chr, b, e, ids]];
  fprintf(' (%i).\n', length(splits)+1);
end

assert(length(unique(chunks(:,4))) == size(chunks,1));
