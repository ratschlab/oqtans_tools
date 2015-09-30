function ev = exon_eval_nostrand(anno_genes, pred_genes)
% ev = intron_eval(anno_genes, pred_genes)
warning('Works only for chromosome size less than 100 mbp!');
warning('UNDER DEVELOPMENT! NOT WORKING CURRENTLY!');

USE_STRAND = 0;

% convert gene annotation and pred annotation into exon lists
anno_exons = extract_unique_exons(anno_genes,USE_STRAND);
pred_exons = extract_unique_exons(pred_genes,USE_STRAND);

anno_exon_matches = zeros(size(anno_exons,1),1);
pred_exon_matches = zeros(size(pred_exons,1),1);

%strands = [+1, -1];
%for s=1:length(strands),
  % TODO hack: this only works for a chromosome size <= 100 Mbp
  %s_idx_anno = find(anno_introns(:,2) == strands(s));
  %anno_intron_ids = anno_introns(s_idx_anno,3:4) + repmat(anno_introns(s_idx_anno,1).*10^8,1,2);
  %s_idx_pred = find(pred_introns(:,2) == strands(s));
  %pred_intron_ids = pred_introns(s_idx_pred,3:4) + repmat(pred_introns(s_idx_pred,1).*10^8,1,2);

  s_idx_anno = 1:size(anno_exons,1);
  s_idx_pred = 1:size(pred_exons,1);
  
  anno_exon_ids = anno_exons(s_idx_anno,3:4) + repmat(anno_exons(s_idx_anno,1).*10^8,1,2);
  pred_exon_ids = pred_exons(s_idx_pred,3:4) + repmat(pred_exons(s_idx_pred,1).*10^8,1,2);

  assert(issorted(anno_exon_ids, 'rows'));
  assert(issorted(pred_exon_ids, 'rows'));

  tol = 0;
  [idx_anno idx_pred] = compare_intervals_sorted(anno_exon_ids, pred_exon_ids, 'similar', tol);
  anno_exon_matches(s_idx_anno) = ~cellfun(@isempty, idx_anno);
  pred_exon_matches(s_idx_pred) = ~cellfun(@isempty, idx_pred);
%end

assert(length(anno_exon_matches) == size(anno_exons,1))
assert(length(pred_exon_matches) == size(pred_exons,1))

fprintf('Compared %i annotated and %i predicted exons\n', size(anno_exons,1), ...
        size(pred_exons,1));
ev.sens = mean(anno_exon_matches);
ev.prec = mean(pred_exon_matches);
ev.f1 = 2 * ev.sens*ev.prec ./ (ev.sens+ev.prec);





function exons = extract_unique_exons(genes,isStrandSpecific)

USE_STRAND = 1;
if exist('isStrandSpecific','var'), 
    USE_STRAND = isStrandSpecific;
end

cnt = 0;
for g=1:length(genes),
  for t=1:length(genes(g).transcripts),
    exons = genes(g).exons{t};
    cnt = cnt + size(exons,1) - 1;
  end
end
exons = zeros(cnt,4);

cnt = 1;
for i=1:length(genes),
    
    % chromosome
    chr = genes(i).chr_num;

    %  discard the strand information?
    strand = genes(i).strand;
    if ~USE_STRAND, strand = '.'; end;

    for t=1:length(genes(i).transcripts),
        % exon list for the current gene transcript
        es = genes(i).exons{t};
        len = size(es,1);

        % add to list
        exons(cnt:cnt+len-1,1) = chr;
        exons(cnt:cnt+len-1,2) = strand;
        exons(cnt:cnt+len-1,3:4) = es;

        cnt = cnt+len;
    end
end

t = size(exons,1);
exons = unique(exons, 'rows');
u = size(exons,1);
fprintf('%i total and %i unique exons\n', t, u);

