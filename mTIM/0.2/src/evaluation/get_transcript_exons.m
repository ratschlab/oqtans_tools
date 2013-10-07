function transcr_exons = get_transcript_exons(genes)
% This file is an adapted version of the get_trancript_introns
% skript by Georg Zeller.
% written by Nico Goernitz, TU Berlin, 2013

cnt = 0;
for g=1:length(genes),
  cnt = cnt + length(genes(g).transcripts);
end
transcr_exons = cell(1,cnt);

cnt = 0;
for g=1:length(genes),
  exons = {};
  for t=1:length(genes(g).transcripts),
    es = genes(g).exons{t};
    num_i = size(es,1);
    s = +1;
    if genes(g).strand == '-', s=-1; end

    exons{end+1} = [repmat([genes(g).chr_num, s], num_i, 1), es];
  end

  if length(exons) >= 1,
    % reduce to the set of alternatively spliced transcripts (i.e. only
    % keep transcripts that differ in their intron structure)
    rm_idx = [];
    for i=2:length(exons),
      for j=1:i-1,
%        fprintf('comparing transcripts %i and %i of gene %i\n', i, j, g);
        % strand and chromosome should always be the same
        assert(unique(exons{i}(:,1)) == unique(exons{j}(:,1)));
        assert(unique(exons{i}(:,2)) == unique(exons{j}(:,2)));
        if isequal(exons{i}, exons{j}),
          rm_idx = [rm_idx j];
        end
      end
    end
%    if ~isempty(rm_idx), 
%      fprintf('removing (intron-wise) redundant transcripts from gene %i\n', g);
%      keyboard
%    end
    exons(unique(rm_idx)) = [];

    transcr_exons(cnt+1:cnt+length(exons)) = exons;
    cnt = cnt + length(exons);
  end
end
assert(cnt <= length(transcr_exons));
transcr_exons(cnt+1:end) = [];

% sort introns
key = zeros(length(transcr_exons),1);
for i=1:length(transcr_exons),
  key(i) = transcr_exons{i}(1,1)*10^8 ...
           + min([transcr_exons{i}(1,3), transcr_exons{i}(end,4)]);
end
[key idx] = sort(key);
transcr_exons = transcr_exons(idx);


