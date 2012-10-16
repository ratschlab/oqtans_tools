function feats = append_mate_feats(chunks, feats, CFG)

% feats = append_mate_feats(chunks, feats, CFG)
%
% 
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009

assert(issorted(chunks, 'rows'));

strands = '+-';
% retrieve position and coverage by interval_query 
% (see interval_query.c in ~/svn/projects/genefinding/utils)
fields = {'Conf_cum'};

L = sum(chunks(:,3)-chunks(:,2)+1);
ssp_feat   = zeros(4, L);

signals = {'don', 'acc'};
for i=1:length(signals),
  data_dir = [CFG.splice_site_dir signals{i} '_pred.bspf/pred/'];
  fprintf('  collecting splice site predictions for %s ...\n', signals{i});
  
  offset = 0;
  tic
  for c=1:size(chunks,1),
    chr = chunks(c,1);
    for s=1:length(strands),
      feat_idx = 2*(s-1)+i;
      fn = sprintf('%scontig_%i%s', data_dir, chr, strands(s));
      [pos, score] = interval_query(fn, fields, [chunks(c,2); chunks(c,3)]);
      switch i,
       case 1,
        % shift donor scores so that they correspond to the first intron position
        switch s,
         case 1,
          pos = pos;
         case 2,
          pos = pos - 2;
        end
       case 2,
        % shift acceptor scores so that they correspond to the first intron position
        switch s,
         case 1,
          pos = pos - 1;
         case 2,
          pos = pos + 1;
        end
      end
      % shifting can result in positions outside the chunk, get rid of those!
      rm_idx = find(pos < chunks(c,2) | chunks(c,3) < pos);
      pos(rm_idx) = [];
      score(rm_idx) = [];
      pos = pos' - chunks(c,2) + 1 + offset;
      ssp_feat(feat_idx, pos) = ssp_feat(feat_idx, pos) + score';
    end
    offset = offset + chunks(c,3)-chunks(c,2)+1;
    if mod(c,100)==0,
      fprintf('    processed %i chunks (%2.1f%%, %.1f sec)\r', ...
              c, 100*c/size(chunks,1), toc);
    end
  end
  fprintf('    processed %i chunks (%2.1f%%, %.1f sec)\n', ...
          c, 100*c/size(chunks,1), toc);
end

assert(~any(any(isnan(ssp_feat))));
feats = [feats; ssp_feat];
assert(size(feats,2) == L);

% eof
