function ev = transcript_eval(anno_genes, pred_genes, use_exons)
% ev = transcript_eval(anno_genes, pred_genes)

USE_EX = 0;
if exist('use_exons','var'),
    USE_EX=use_exons; 
end

if ~USE_EX,
    anno_transcr = get_transcript_introns(anno_genes);
    pred_transcr = get_transcript_introns(pred_genes);
else
    anno_transcr = get_transcript_exons(anno_genes);
    pred_transcr = get_transcript_exons(pred_genes);
end

% chose this to evaluate only a single randomly chosen transcript per
% predicted gene (discarding all alternative transcripts)
%pred_transcr = get_transcript_introns_no_alt(pred_genes);

anno_transcr_strand = zeros(length(anno_transcr),1);
for i=1:length(anno_transcr),
  anno_transcr_strand(i) = anno_transcr{i}(1,2);
end
pred_transcr_strand = zeros(length(pred_transcr),1);
for i=1:length(pred_transcr),
  pred_transcr_strand(i) = pred_transcr{i}(1,2);
end

anno_transcr_matches = zeros(length(anno_transcr),1);
pred_transcr_matches = zeros(length(pred_transcr),1);

strands = [+1, -1];
str = {'+', '-'};
for s=1:length(strands),
  s_idx_anno = find(anno_transcr_strand == strands(s));
  s_idx_pred = find(pred_transcr_strand == strands(s));

  % First, determine transcript boundaries and overlapping transcripts
  % TODO hack: this only works for a chromosome size <= 100 Mbp
  anno_boundary_ids = zeros(length(s_idx_anno),2);
  for i=1:length(s_idx_anno),
    t = s_idx_anno(i);
    anno_boundary_ids(i,1) = anno_transcr{t}(1,3);
    assert(all(anno_boundary_ids(i,1) <= anno_transcr{t}(:,3)));
    anno_boundary_ids(i,1) = anno_boundary_ids(i,1) + anno_transcr{t}(1,1) * 10^8;

    anno_boundary_ids(i,2) = anno_transcr{t}(end,4);
    assert(all(anno_boundary_ids(i,2) >= anno_transcr{t}(:,4)));
    anno_boundary_ids(i,2) = anno_boundary_ids(i,2) + anno_transcr{t}(1,1) * 10^8;
  end
  pred_boundary_ids = zeros(length(s_idx_pred),2);
  for i=1:length(s_idx_pred),
    t = s_idx_pred(i);
    pred_boundary_ids(i,1) = pred_transcr{t}(1,3);
    assert(all(pred_boundary_ids(i,1) <= pred_transcr{t}(:,3)));
    pred_boundary_ids(i,1) = pred_boundary_ids(i,1) + pred_transcr{t}(1,1) * 10^8;

    pred_boundary_ids(i,2) = pred_transcr{t}(end,4);
    assert(all(pred_boundary_ids(i,2) >= pred_transcr{t}(:,4)));
    pred_boundary_ids(i,2) = pred_boundary_ids(i,2) + pred_transcr{t}(1,1) * 10^8;
  end

%  tic
  [anno_boundary_ids anno_perm] = sortrows(anno_boundary_ids);
  s_idx_anno = s_idx_anno(anno_perm);
  [pred_boundary_ids pred_perm] = sortrows(pred_boundary_ids);
  s_idx_pred = s_idx_pred(pred_perm);
%  assert(issorted(anno_boundary_ids, 'rows'));
%  assert(issorted(pred_boundary_ids, 'rows'));
  ovl = 1;
  [ovl_anno ovl_pred] = compare_intervals_sorted(anno_boundary_ids, pred_boundary_ids, ...
                                                 'overlap', ovl);
%  fprintf(['Interval comparison between annotated and predicted transcript ' ...
%           'boundaries took %.1f sec\n'], toc);
  % Second, compare the intron structure of all of their transcripts
  tic
  for i=1:length(ovl_anno),
    match = 0;
    for j=1:length(ovl_anno{i}),
      k = ovl_anno{i}(j);
      if isequal(anno_transcr{s_idx_anno(i)}, pred_transcr{s_idx_pred(k)})
        match = 1;
        break
      end
    end
    anno_transcr_matches(s_idx_anno(i)) = match;
  end
%  fprintf('Sensitivity assessment took %.1f sec                          \n', toc);
%  fprintf('Mean sensitivity on %s strand: %2.2f%%\n', str{s}, 100* ...
%          mean(anno_transcr_matches(s_idx_anno)));

  tic
  for i=1:length(ovl_pred),
    match = 0;
    for j=1:length(ovl_pred{i})
      k = ovl_pred{i}(j);
      if isequal(pred_transcr{s_idx_pred(i)}, anno_transcr{s_idx_anno(k)})
        match = 1;
        break
      end
    end
    pred_transcr_matches(s_idx_pred(i)) = match;
  end
%  fprintf('Precision assessment took %.1f sec                            \n', toc);
%  fprintf('Mean precision on %s strand: %2.2f%%\n', str{s}, 100* ...
%          mean(pred_transcr_matches(s_idx_pred)));

end
fprintf('Compared %i annotated and %i predicted transcripts\n', ...
        length(anno_transcr), length(pred_transcr));

ev.sens = mean(anno_transcr_matches);
ev.prec = mean(pred_transcr_matches);
ev.f1 = 2 * ev.sens*ev.prec ./ (ev.sens+ev.prec);

