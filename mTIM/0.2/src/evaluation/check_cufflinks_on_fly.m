function check_cufflinks_on_fly(CFG)
% IN:
%   CFG     : (structure) containing directories and files
%             for evaluation and corresponding label
%
warning('This is NOT a valid function! Only for the paper!');

anno_dir = CFG.annotation_dir;

% Open text output file:
% 1. if file doesn't exist create a new one
% 2. if file already exist discard the content 
fh = fopen('/home/nico/0/pred.txt','w');


% load cufflinks (or other competitors) predicted gene structure
cuffl = [];
CFG.cufflinks_pred_file = '/home/nico/0/cuffl.mat';
fprintf('Loading cufflinks genes from "%s"...\n', CFG.cufflinks_pred_file);
fprintf(fh,'Loading cufflinks genes from "%s"...\n', CFG.cufflinks_pred_file);

cuffl = load(CFG.cufflinks_pred_file, 'genes');

fprintf('  loaded %i predicted genes.\n', length(cuffl.genes));
fprintf(fh,'  loaded %i predicted genes.\n', length(cuffl.genes));
fprintf('  converting genes: %i\n',CFG.cufflinks_convert);
cuffl.genes = closed_to_half_open(cuffl.genes);
fprintf('Done!\n');

% load annotation
fprintf('Loading annotation...\n');
anno = load([anno_dir 'genes.mat'], 'genes');
fprintf('  using %i annotated genes for evaluation.\n', length(anno.genes));
fprintf(fh,'  using %i annotated genes for evaluation.\n', length(anno.genes));
anno.genes = prune_gene_struct(anno.genes);


%%% Sort gene structures (prerequisite for fast comparisons on sorted
%intervals in the following)
anno.genes = merge_genes_by_name_elegans(anno.genes);
gid = 10^8*[anno.genes.chr_num] + min([[anno.genes.start]; [anno.genes.stop]]);
[tmp perm] = sort(gid);
anno.genes = anno.genes(perm);

%keyboard
%new.genes = split_genes_filter(CFG,pred.genes);
for i=1:length(cuffl.genes),
    new.genes(i) = cuffl.genes(i);
    nt = length(cuffl.genes(i).transcripts);
    ind = randperm(nt);
    new.genes(i).transcripts = cuffl.genes(i).transcripts(ind(1));
end
[eval, split, fusion] = evaluate(anno.genes,new.genes,fh);


