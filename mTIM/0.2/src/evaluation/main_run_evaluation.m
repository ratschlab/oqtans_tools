function main_run_evaluation(CFG)
% IN:
%   CFG     : (structure) containing directories and files
%             for evaluation and corresponding label
%
% 
%

anno_dir = CFG.annotation_dir;
pred_dir = CFG.model_dirs{CFG.opt_model};
fname = CFG.pred_filename;
out_prefix = 'eval_prediction';
if (isfield(CFG,'eval_vald') && CFG.eval_vald),
    out_prefix = 'eval_validation';
    fname = CFG.vald_filename;
    fprintf('Evaluate validation data.\n');
end

% save evaluation into specific file
if isfield(CFG,'eval_filename'),
    out_prefix = CFG.eval_filename;
    fprintf('Use %s as evaluation output file name.\n');
end


% Open text output file:
% 1. if file doesn't exist create a new one
% 2. if file already exist discard the content 
fh = fopen([pred_dir out_prefix '.txt'],'w');


% load cufflinks (or other competitors) predicted gene structure
cuffl = [];
if (isfield(CFG,'cufflinks_pred_file')),
    fprintf('Loading cufflinks genes from "%s"...\n', CFG.cufflinks_pred_file);
    fprintf(fh,'Loading cufflinks genes from "%s"...\n', CFG.cufflinks_pred_file);
    cuffl = load(CFG.cufflinks_pred_file, 'genes');
    fprintf('  loaded %i predicted genes.\n', length(cuffl.genes));
    fprintf(fh,'  loaded %i predicted genes.\n', length(cuffl.genes));
    fprintf('  converting genes: %i\n',CFG.cufflinks_convert);

    if (CFG.cufflinks_convert),
        cuffl.genes = closed_to_half_open(cuffl.genes);
    end
    fprintf('Done!\n');
end

% load predictions
fprintf('Loading predicted genes (%s).\n', [pred_dir fname]);
fprintf(fh,'Loading predicted genes (%s).\n', [pred_dir fname]);
pred = load([pred_dir fname], 'genes');
fprintf('  evaluating %i predicted genes.\n', length(pred.genes));
fprintf(fh,'  evaluating %i predicted genes.\n', length(pred.genes));


% load annotation
fprintf('Loading annotation...\n');
anno = load([anno_dir 'genes.mat'], 'genes');
fprintf('  using %i annotated genes for evaluation.\n', length(anno.genes));
fprintf(fh,'  using %i annotated genes for evaluation.\n', length(anno.genes));
anno.genes = prune_gene_struct(anno.genes);

% convert genes
pred.genes = closed_to_half_open(pred.genes);

% remove all genes with too low coverage
%idxs = find([pred.genes.expr]>1.1);
%fprintf('Remove %i of %i genes due to low coverage.\n',length(pred.genes)-length(idxs), length(pred.genes));
%pred.genes = pred.genes(idxs);



%%% Sort gene structures (prerequisite for fast comparisons on sorted
%intervals in the following)
anno.genes = merge_genes_by_name_elegans(anno.genes);
gid = 10^8*[anno.genes.chr_num] + min([[anno.genes.start]; [anno.genes.stop]]);
[tmp perm] = sort(gid);
anno.genes = anno.genes(perm);
gid = 10^8*[pred.genes.chr_num] + min([[pred.genes.start]; [pred.genes.stop]]);
[tmp perm] = sort(gid);
pred.genes = pred.genes(perm);

%keyboard
%new.genes = split_genes_filter(CFG,pred.genes);
%for i=1:length(anno.genes),
%    new.genes(i) = anno.genes(i);
%    new.genes(i).transcripts = anno.genes(i).transcripts(1);
%end
%[eval, split, fusion] = evaluate(anno.genes,new.genes,fh);

% mTIM
fprintf('\n mTIM results:\n');
fprintf(fh, '\n mTIM results:\n');
[eval, split, fusion] = evaluate(anno.genes, pred.genes, fh);

fprintf('Filter splitted and merged genes from data set:\n');
fprintf(fh,'Filter splitted and merged genes from data set:\n');

anno_inds = setdiff(1:length(anno.genes),[split.anno_inds, fusion.anno_inds]);
pred_inds = setdiff(1:length(pred.genes),[split.pred_inds, fusion.pred_inds]);
[eval_filter, split, fusion] = evaluate(anno.genes(anno_inds), pred.genes(pred_inds), fh);

% CUFFLINKS
fprintf('\n Cufflinks results:\n');
fprintf(fh, '\n Cufflinks results:\n');

if (~isempty(cuffl)),
    [eval_cuffl, split, fusion] = evaluate(anno.genes, cuffl.genes, fh);

    fprintf('Filter splitted and merged genes from data set:\n');
    fprintf(fh,'Filter splitted and merged genes from data set:\n');

    % vanilla genes (no splits and no fusion)
    anno_inds = setdiff(1:length(anno.genes),[split.anno_inds, fusion.anno_inds]);
    cuffl_inds = setdiff(1:length(cuffl.genes),[split.pred_inds, fusion.pred_inds]);
    [eval_cuffl_filter, split, fusion] = evaluate(anno.genes(anno_inds), cuffl.genes(cuffl_inds), fh);

    % save result
    save([pred_dir out_prefix '.mat'], 'eval','eval_cuffl');
else
    fprintf('\n No Cufflinks files defined.\n');
    fprintf(fh, '\n No Cufflinks files defined.\n');
    save([pred_dir out_prefix '.mat'], 'eval','eval_filter');
end

% add some information
txt = print_structure(CFG);

fprintf(fh,'\n\nConfiguration:\n%s',txt);
fclose(fh);
