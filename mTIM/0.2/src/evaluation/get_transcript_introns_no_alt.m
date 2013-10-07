function transcr_introns = get_transcript_introns_no_alt(genes)

cnt = 0;
for g=1:length(genes),
  cnt = cnt + 1;
end
transcr_introns = {};

for g=1:length(genes),
  % choose one transcript at random
  t = ceil(length(genes(g).transcripts).*rand(1));
  exons = genes(g).exons{t};
  num_i = size(exons,1) - 1;
  s = +1;
  if genes(g).strand == '-',
    s=-1;
  end
  introns = [repmat([genes(g).chr_num, s], num_i, 1), ...
             exons(1:end-1,2), exons(2:end,1)];
  % do not consider single-exon transcripts
  if ~isempty(introns),
    transcr_introns{end+1} = introns;
  end
end

% sort introns
key = zeros(length(transcr_introns),1);
for i=1:length(transcr_introns),
  key(i) = transcr_introns{i}(1,1)*10^8 ...
           + min([transcr_introns{i}(1,3), transcr_introns{i}(end,4)]);
end
[key idx] = sort(key);
transcr_introns = transcr_introns(idx);


