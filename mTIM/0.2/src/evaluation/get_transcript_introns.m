function transcr_introns = get_transcript_introns(genes)

cnt = 0;
for g=1:length(genes),
  cnt = cnt + length(genes(g).transcripts);
end
transcr_introns = cell(1,cnt);

cnt = 0;
for g=1:length(genes),
  introns = {};
  for t=1:length(genes(g).transcripts),
    exons = genes(g).exons{t};
    num_i = size(exons,1) - 1;
    s = +1;
    if genes(g).strand == '-',
      s=-1;
    end
    introns{end+1} = [repmat([genes(g).chr_num, s], num_i, 1), ...
                      exons(1:end-1,2), exons(2:end,1)];
  end
  % remove single-exon transcripts
  rm_idx = cellfun(@isempty, introns);
  introns(rm_idx) = [];

  if length(introns) >= 1,
    % reduce to the set of alternatively spliced transcripts (i.e. only
    % keep transcripts that differ in their intron structure)
    rm_idx = [];
    for i=2:length(introns),
      for j=1:i-1,
%        fprintf('comparing transcripts %i and %i of gene %i\n', i, j, g);
        % strand and chromosome should always be the same
        assert(unique(introns{i}(:,1)) == unique(introns{j}(:,1)));
        assert(unique(introns{i}(:,2)) == unique(introns{j}(:,2)));
        if isequal(introns{i}, introns{j}),
          rm_idx = [rm_idx j];
        end
      end
    end
%    if ~isempty(rm_idx), 
%      fprintf('removing (intron-wise) redundant transcripts from gene %i\n', g);
%      keyboard
%    end
    introns(unique(rm_idx)) = [];

    transcr_introns(cnt+1:cnt+length(introns)) = introns;
    cnt = cnt + length(introns);
  end
end
assert(cnt <= length(transcr_introns));
transcr_introns(cnt+1:end) = [];

% sort introns
key = zeros(length(transcr_introns),1);
for i=1:length(transcr_introns),
  key(i) = transcr_introns{i}(1,1)*10^8 ...
           + min([transcr_introns{i}(1,3), transcr_introns{i}(end,4)]);
end
[key idx] = sort(key);
transcr_introns = transcr_introns(idx);


