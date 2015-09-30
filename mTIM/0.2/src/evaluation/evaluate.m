function [eval, split, fusion] = evaluate(ag, pg, fh)
% ag - annotation genes
% pg - prediction genes
% fh - (optional) file handle to write text output to

%%% Intron-level evaluation
eval = struct();
eval.exon = exon_eval_nostrand(ag, pg);
fprintf('Exon evaluation (w/o strand):\n');
eval.exon

%%% Intron-level evaluation
eval = struct();
eval.intron = intron_eval(ag, pg);
fprintf('Intron evaluation:\n');
eval.intron

%%% Intron-based transcript-level evaluation
eval.transcr = transcript_eval(ag, pg);
fprintf('Transcript evaluation:\n');
eval.transcr


%%% Intron-based gene-level evaluation
[eval.gene, split, fusion] = gene_eval(ag, pg);
fprintf('Gene evaluation:\n');
eval.gene


fprintf('SPLITTING: There are %i in annotation and %i genes in prediction effected.\n', ...
    length(split.anno_inds),length(split.pred_inds));
fprintf('MERGING: There are %i in annotation and %i genes in prediction effected.\n', ...
    length(fusion.anno_inds),length(fusion.pred_inds));


% if a file handle is present then write
% evaluation output in text format
if (exist('fh','var')),
    fprintf(fh,'SPLITTING: There are %i in annotation and %i genes in prediction effected.\n', ...
        length(split.anno_inds),length(split.pred_inds));
    fprintf(fh,'MERGING: There are %i in annotation and %i genes in prediction effected.\n', ...
        length(fusion.anno_inds),length(fusion.pred_inds));

    fprintf(fh, 'Intron evaluation:\n');
    fprintf(fh, ' sens: %1.4f\n',eval.intron.sens);
    fprintf(fh, ' prec: %1.4f\n',eval.intron.prec);
    fprintf(fh, ' f1  : %1.4f\n',eval.intron.f1);

    fprintf(fh, 'Transcript evaluation:\n');
    fprintf(fh, ' sens: %1.4f\n',eval.transcr.sens);
    fprintf(fh, ' prec: %1.4f\n',eval.transcr.prec);
    fprintf(fh, ' f1  : %1.4f\n',eval.transcr.f1);

    fprintf(fh, 'Gene evaluation:\n');
    fprintf(fh, ' sens: %1.4f\n',eval.gene.sens);
    fprintf(fh, ' prec: %1.4f\n',eval.gene.prec);
    fprintf(fh, ' f1  : %1.4f\n',eval.gene.f1);
end

