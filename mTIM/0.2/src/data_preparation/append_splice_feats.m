function feats = append_splice_feats(chunks, feats, CFG)
% feats = append_splice_feats(chunks, feats, CFG)
%
% 
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2009-2011
%            Nico Goernitz, TU Berlin, 2012 

if ~isfield(CFG.PAR,'use_only_consensus_splice_sites')
    % *DEFAULT*
    % if there is no mention of splice feature type at 
    % all then use computational splice sites
    feats = append_comp_splice_sites(chunks, feats, CFG);
else

    if CFG.PAR.use_only_consensus_splice_sites,
        feats = append_splice_consensus_feats(chunks,feats,CFG);
    else
        feats = append_comp_splice_sites(chunks, feats, CFG);
    end

end



function feats = append_splice_consensus_feats(chunks, feats, CFG)
% Mark only consensus (if no computational splice sites are available).
%
%
warning('SPLICE CONSENSUS FEATS ARE NOT TESTED YET.');
warning('WITHOUT ANY CONSENSUS.');

assert(issorted(chunks, 'rows'));
tic

L = sum(chunks(:,3)-chunks(:,2)+1);
%ssp_feat = -inf(4, L);
ssp_feat = zeros(4, L);


genome_info = init_genome(CFG.genome_info);
for c=1:CFG.num_chr
  genome_info.flat_fnames{c} = [genome_info.basedir '/genome/' genome_info.contig_names{c} '.flat'];
  seq{c} = load_genomic(c,'+',1,CFG.chr_lens(c), genome_info);
end

pos = 0;
for c=1:size(chunks,1),
    chr = chunks(c,1);
    s_start = chunks(c,2);
    s_end = chunks(c,3);
    s = seq{chr}(s_start:s_end);
    % search for consensus strings within 's'
    % W-strand
    %   donor: 'gt' 'gc'
    %   acceptor: 'ag'
    % C-strand
    %   donor: 'ct'
    %   acceptor: 'ac'

    % format [don_W acc_W don_C acc_C]
    ssp_feat(1,pos + findstr(s,'gt')) = 1;
    ssp_feat(2,pos + findstr(s,'ag')) = 1;
   
    ssp_feat(3,pos + findstr(s,'tg')) = 1;
    ssp_feat(4,pos + findstr(s,'ga')) = 1;
  
    pos = pos + s_end-s_start + 1;
end

assert(~any(any(isnan(ssp_feat))));
feats = [feats; ssp_feat];
assert(size(feats,2) == L);

fprintf('  appended splice site CONSENSUS features to %i chunks in %.1f sec\n', ...
        size(chunks,1), toc);







function feats = append_comp_splice_sites(chunks,feats,CFG)
% Append computational splice sites as features if 
% available.
%
assert(issorted(chunks, 'rows'));
tic

strands = '+-';
% retrieve position and coverage by interval_query 
% (see interval_query.c in ~/svn/projects/genefinding/utils)
fields = {'Conf_cum'};

L = sum(chunks(:,3)-chunks(:,2)+1);
ssp_feat = -inf(4, L);

signals = {'don', 'acc'};
for i=1:length(signals),
  data_dir = [getfield(CFG.splice_site_dir, signals{i}) '/'];
  if exist([data_dir 'pred/'], 'dir'),
    data_dir = [getfield(CFG.splice_site_dir, signals{i}) '/pred/'];
  end
  if isempty([data_dir '*.pos'])
    error('Failed to locate splice site predictions at %s\n', data_dir);
  end
%  fprintf('Loading splice site predictions from %s\n', data_dir);
%  fprintf('  collecting splice site predictions for %s ...\n', signals{i});
  
  offset = 0;
  for c=1:size(chunks,1),
    chr = chunks(c,1);
    for s=1:length(strands),
      feat_idx = 2*(s-1)+i;
      fn = sprintf('%scontig_%i%s', data_dir, chr, strands(s));
      if ~exist([fn '.pos'], 'file'),
        error('Failed to locate splice site predictions for contig_%i%s at %s\n', ...
              chr, strands(s), data_dir);
      end
      [pos, score] = interval_query(fn, fields, [chunks(c,2); chunks(c,3)]);
      switch i,
       case 1,
        % no shifting needed for donor scores to be at the last intron position
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
      ssp_feat(feat_idx, pos) = score';
    end
    offset = offset + chunks(c,3)-chunks(c,2)+1;
  end
end

assert(~any(any(isnan(ssp_feat))));
feats = [feats; ssp_feat];
assert(size(feats,2) == L);

fprintf('  appended predicted splice site features to %i chunks in %.1f sec\n', ...
        size(chunks,1), toc);

% eof
