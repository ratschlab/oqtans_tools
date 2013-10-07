function main_run_visualization(CFG)
% IN:
%   CFG     : (structure) containing directories and files
%             for evaluation and corresponding label
%
%

VIS_BROWSER = 1;
VIS_BROWSER_XCORR = 0;

VIS_EXPR_METHOD = 0;
VIS_EXPR_BINS = 5;
VIS_EXPR = 0;


anno_dir = CFG.annotation_dir;
pred_dir = CFG.model_dirs{CFG.opt_model};
fname = CFG.pred_filename;


% load cufflinks (or other competitors) predicted gene structure
cuffl = [];
if (isfield(CFG,'cufflinks_pred_file')),
    fprintf('Loading cufflinks genes from "%s"...\n', CFG.cufflinks_pred_file);
    cuffl = load(CFG.cufflinks_pred_file, 'genes');
    fprintf('  loaded %i predicted genes.\n', length(cuffl.genes));
    fprintf('  converting genes: %i\n',CFG.cufflinks_convert);
    if (CFG.cufflinks_convert),
        cuffl.genes = closed_to_half_open(cuffl.genes);
    end
    fprintf('Done!\n');
end
 
% load predictions
fprintf('Loading predicted genes (%s).\n', [pred_dir fname]);
pred = load([pred_dir fname], 'genes');
fprintf('  evaluating %i predicted genes.\n', length(pred.genes));

% load annotation
fprintf('Loading annotation...\n');
anno = load([anno_dir 'genes.mat'], 'genes');
fprintf('  using %i annotated genes for evaluation.\n', length(anno.genes));
anno.genes = prune_gene_struct(anno.genes);

% convert genes
pred.genes = closed_to_half_open(pred.genes);

%%% Sort gene structures (prerequisite for fast comparisons on sorted
%intervals in the following)
anno.genes = merge_genes_by_name_elegans(anno.genes);
gid = 10^8*[anno.genes.chr_num] + min([[anno.genes.start]; [anno.genes.stop]]);
[tmp perm] = sort(gid);
anno.genes = anno.genes(perm);
gid = 10^8*[pred.genes.chr_num] + min([[pred.genes.start]; [pred.genes.stop]]);
[tmp perm] = sort(gid);
pred.genes = pred.genes(perm);

data_fn = [CFG.xval_dirs{1} 'test_data.mat'];
load(data_fn, '-mat', 'test_chunks');
fprintf('Using test chunks from %s\n', data_fn);

%keyboard
%new.genes = split_genes_filter(CFG,pred.genes);
%[eval, split, fusion] = evaluate(anno.genes,new.genes);


%fprintf('\n\nOverall evaluation.\n');
%%% Intron-level evaluation
%eval = struct();
%eval.intron = intron_eval(anno.genes, pred.genes);
%fprintf('Intron evaluation:\n');
%eval.intron

%%% Intron-based transcript-level evaluation
%eval.transcr = transcript_eval(anno.genes, pred.genes);
%fprintf('Transcript evaluation:\n');
%eval.transcr

%%% Intron-based gene-level evaluation
[eval.gene, split, fusion] = gene_eval(anno.genes, pred.genes);
fprintf('Gene evaluation:\n');
eval.gene


%fprintf('Filter splitted and merged genes from data set.\n');

% vanilla genes (no splits and no fusion)
%anno_inds = setdiff(1:length(anno.genes),[split.anno_inds, fusion.anno_inds]);
%pred_inds = setdiff(1:length(pred.genes),[split.pred_inds, fusion.pred_inds]);

%vanilla_inds = intersect(sinds,finds);
%fprintf('SPLITTING: There are %i in annotation and %i genes in prediction effected.\n', ...
%    length(split.anno_inds),length(split.pred_inds));
%fprintf('MERGING: There are %i in annotation and %i genes in prediction effected.\n', ...
%    length(fusion.anno_inds),length(fusion.pred_inds));


% i want to see only genes that have been combined
fprintf('%i predicted genes are affected by fusion.\n',length(fusion.pred_inds));
my_chunks = [];
for i=1:length(fusion.pred_inds),
    ind = fusion.pred_inds(i);

%    if ~isempty(strfind(names{i},'combined')),
%        start = new.genes(i).start - 1000;
%        stop = new.genes(i).stop + 1000;
%        chr = new.genes(i).chr_num;
%        my_chunks = [my_chunks; chr,start,stop,-1,-1];
%    end
     start = pred.genes(ind).start - 1000;
     stop = pred.genes(ind).stop + 1000;
     chr = pred.genes(ind).chr_num;
     my_chunks = [my_chunks; chr,start,stop,-1,-1];

end

if (VIS_BROWSER),
  %show_gene_browser(CFG, test_chunks, anno.genes, pred.genes, new.genes, VIS_BROWSER_XCORR);
  show_gene_browser(CFG, my_chunks, anno.genes, pred.genes, cuffl.genes, VIS_BROWSER_XCORR);
end


if (VIS_EXPR),
    if (~isempty(cuffl)),
      fprintf('Cufflinks:\n');
      levels = get_expression_bins(cuffl.genes,VIS_EXPR_BINS,VIS_EXPR_METHOD);
      plot_levels(levels, anno, cuffl, 'r');
    end
    fprintf('mTIM:\n');
    levels = get_expression_bins(pred.genes,VIS_EXPR_BINS,VIS_EXPR_METHOD);
    plot_levels(levels, anno, pred, 'b');
end




function plot_levels(levels, anno, pred, color);
num_levels = length(levels);
figure(2);
hold on;
for i=1:num_levels,
        inds = levels{i};
        fprintf('There are %i genes within level %i.\n',length(inds),i);

        subplot(2,2,1);
        hold on;
        alpha(0.5);
        bar(i,length(inds));
        hold off;

        % evaluate
        eval = evaluate(anno.genes, pred.genes(inds));

        subplot(2,2,2);
        hold on;
        alpha(0.5);
        bar(i,eval.intron.prec,color);
        hold off;

        subplot(2,2,3);
        hold on;
        alpha(0.5);
        bar(i,eval.transcr.prec,color);
        hold off;

        subplot(2,2,4);
        hold on;
        alpha(0.5);
        bar(i,eval.gene.prec,color);
        hold off;
     
end
hold off;


% eof
