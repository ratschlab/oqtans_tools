function eval_gene_pred()

VIZ_EVAL = 1
LEVEL_EVAL = 0

% TODO hacky
addpath ../utils
addpath ../utils2

%method = 'cufflinks'
%method = 'scripture'
method = 'mtim'

organism = 'elegans'
%organism = 'melanogaster'
%organism = 'thaliana'

%exp = '6';
exp = '6.bestpaired.filtered';
%exp = 'tophat.PE';
%exp = 'tophat.PE.bestpaired.filtered';
%exp = 'SRX_all_new';

switch organism,
  case 'elegans',
   genome_dir = ['/fml/ag-raetsch/nobackup/projects/rgasp.2/genomes/elegans/' ...
                 'elegans.gio/'];
   
   anno_dir = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/elegans/';

   data_dir = '/fml/ag-raetsch/nobackup/projects/sequencing_runs/worm-Dec09/reads/mTIM/';
   switch method,
    case 'cufflinks',
     pred_dir = sprintf('%scufflinks/elegans/elegans.%s/sanitized_out/', data_dir, exp);
    case 'scripture',
     pred_dir = sprintf('%sscripture/elegans/elegans.%s/', data_dir, exp);
    case 'mtim',
%     exp_dir = 'mtim_pair_10Mar'; 
%     pred_dir = sprintf('%s%s.elegans.%s/galaxy/predictions/', data_dir, exp_dir, exp);
     
     % the prediction .mat-file is within this directory
     pred_dir = '';
   end
 case 'melanogaster',
   genome_dir = ['/fml/ag-raetsch/nobackup/projects/rgasp.2/genomes/drosophila/' ...
                 'drosophila.gio/'];
   
   anno_dir = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/drosophila/';

   data_dir = '/fml/ag-raetsch/nobackup/projects/sequencing_runs/D_melanogaster/reads/';
   switch method,
    case 'cufflinks',
     pred_dir = sprintf('%scufflinks/drosophila/D_melanogaster_L3.%s/sanitized_out/', data_dir, exp);
    case 'scripture',
     pred_dir = sprintf('%sscripture/drosophila/D_melanogaster_L3.%s/', data_dir, exp);
    case 'mtim',
%     exp_dir = 'mtim_exp_pair_04Apr';
%     exp_dir = 'mtim_exp_pair_02May';
     exp_dir = 'mtim_exp_pair_11May';
%     exp_dir = 'mtim_exp_pair_12May';
     pred_dir = sprintf('%s%s.melanogaster.%s/galaxy/predictions/', data_dir, exp_dir, exp);
   end
 case 'thaliana',
   data_dir = '/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana/';
   genome_dir = [data_dir 'genome_info_for_mTIM/'];
   anno_dir = [data_dir 'annotation_for_mTIM/'];
  switch method,
   case 'cufflinks',
    pred_dir = [data_dir 'cufflinks/'];
   case 'scripture',
    error('no scripture predictions available for A. thaliana');
   case 'mtim',
%     exp_dir = 'mtim_pair_12Mar'; 
     pred_dir = sprintf('%s%s.%s/galaxy/predictions/', data_dir, exp_dir, exp);
  end
 case 'musculus',
  error('not yet implemented');
end

switch method,
 case 'mtim',
% C. elegans
%

% D. melanogaster
%  train_time = '2011-04-01_11h54', pred_time = '2011-04-04_09h25' % 12.Mar.melanogaster.6.bestpaired.filtered 
%  train_time = '2011-04-01_11h54', pred_time = '2011-04-04_10h01' % 12.Mar.melanogaster.6.bestpaired.filtered 
%  train_time = '2011-04-12_09h01', pred_time = '2011-04-13_10h36' % 04Apr.melanogaster.6.bestpaired.filtered 
%  train_time = '2011-04-12_13h31', pred_time = '2011-04-13_10h39' % 04Apr.melanogaster.6.bestpaired.filtered 
%  train_time = '2011-04-13_13h10', pred_time = '2011-04-28_14h56' % 04Apr.melanogaster.6.bestpaired.filtered 
%  train_time = '2011-04-13_13h10', pred_time = '2011-05-09_14h02' % 04Apr.melanogaster.6.bestpaired.filtered after gen_graph
%  train_time = '2011-05-03_08h07', pred_time = '2011-05-04_07h37' % 02May.melanogaster.6.bestpaired.filtered 
%  train_time = '2011-05-03_08h42', pred_time = '2011-05-04_07h38' % 02May.melanogaster.6.bestpaired.filtered 
%  train_time = '2011-05-03_09h13', pred_time = '2011-05-06_13h53' % 02May.melanogaster.6.bestpaired.filtered 
%  train_time = '2011-05-11_11h00', pred_time = '2011-05-12_08h56' % 11May.melanogaster.6.bestpaired.filtered
%  train_time = '2011-05-11_12h08', pred_time = '2011-05-12_08h50' % 11May.melanogaster.6.bestpaired.filtered
  train_time = '2011-05-12_14h02', pred_time = '2011-05-16_15h17' % 11May.melanogaster.6.bestpaired.filtered
%  train_time = '2011-05-12_14h18', pred_time = '2011-05-16_14h41' % 12May.melanogaster.6.bestpaired.filtered

% A. thaliana
%

  fname = sprintf('genes_pred%s_trained%s', pred_time, train_time);
  warning('file name hack.');
  fname = 'pred'
 case 'cufflinks',
  fname = 'transcripts.gtf';
 case 'scripture',
  fname = 'sanitized.segments';  
 otherwise
  error('unknown method: %s', method);
end

% load predictions
pred_dir = '~/Projects/new_mTIM/out/elegans/120125_1534/model1/';
fname = 'prediction';
fprintf('Loading predicted genes (%s).\n', [pred_dir fname]);
pred = load([pred_dir fname '.mat'], 'genes');
fprintf('  evaluating %i predicted genes.\n', length(pred.genes));

% initialize genome config
fn_genome_config = [genome_dir 'genome.config'];
genome_info = init_genome(fn_genome_config);

% load annotation
fprintf('Loading annotation...\n');
anno = load([anno_dir 'genes.mat'], 'genes');
fprintf('  using %i annotated genes for evaluation.\n', length(anno.genes));
anno.genes = prune_gene_struct(anno.genes);

switch method,
 case 'mtim',
  % convert gene structure from closed intervals to half-open intervals
% TODO
  pred.genes = closed_to_half_open(pred.genes);
 case 'cufflinks',
  % convert gene structure from closed intervals to half-open intervals
  pred.genes = closed_to_half_open(pred.genes);
 case 'scripture'
  % no need to convert these
  ;
end


%keyboard
% overwrite pred with RGASP.2 submission here if desired
% C. elegans
%pred = load(['/fml/ag-raetsch/nobackup/projects/rgasp.2/submissions/' ...
%             'elegans/elegans_C_elegans_L3_mTIM:filtered_gff.mat'], 'genes');
% D. melanogaster
%pred = load(['/fml/ag-raetsch/nobackup/projects/rgasp.2/submissions/' ...
%             'drosophila/drosophila_D_melanogaster_L3_mTIM:filtered_gff.mat'], 'genes');


%keyboard
% remove alternative transcripts randomly retaining one isoform
%for g=1:length(pred.genes),
%  r = ceil(length(pred.genes(g).transcripts).*rand(1));
%  idx = setdiff(1:length(pred.genes(g).transcripts), r);
%  pred.genes(g).transcripts(idx) = [];
%  pred.genes(g).exons(idx) = [];
%  pred.genes(g).start = pred.genes(g).exons{1}(1,1);
%  pred.genes(g).stop = pred.genes(g).exons{1}(end,2);
%end

%%% Sort gene structures (prerequisite for fast comparisons on sorted
%intervals in the following)
anno.genes = merge_genes_by_name_elegans(anno.genes);
gid = 10^8*[anno.genes.chr_num] + min([[anno.genes.start]; [anno.genes.stop]]);
[tmp perm] = sort(gid);
anno.genes = anno.genes(perm);
gid = 10^8*[pred.genes.chr_num] + min([[pred.genes.start]; [pred.genes.stop]]);
[tmp perm] = sort(gid);
pred.genes = pred.genes(perm);

%%% Intron-level evaluation
eval = struct();
eval.intron = intron_eval(anno.genes, pred.genes);
fprintf('Intron evaluation:\n');
eval.intron

%%% Intron-based transcript-level evaluation
eval.transcr = transcript_eval(anno.genes, pred.genes);
fprintf('Transcript evaluation:\n');
eval.transcr

%%% Intron-based gene-level evaluation
eval.gene = gene_eval(anno.genes, pred.genes);
fprintf('Gene evaluation:\n');
eval.gene

eval_dir = sprintf('%smtim_eval/%s/%s/', data_dir, exp, method);
if ~exist(sprintf('%smtim_eval/%s', data_dir, exp)),
  mkdir(sprintf('%smtim_eval/%s', data_dir, exp));
end
if ~exist(eval_dir),
  mkdir(eval_dir)
end
ev_fname = ['eval_' fname];
fprintf('HACK.. save file to the current directory.\n');
%save([eval_dir ev_fname '.mat'], 'eval');
save([ev_fname '.mat'], 'eval');

%%% Visualization of predictions
if isequal(method, 'mtim') && VIZ_EVAL,
  switch organism,
   case 'melanogaster',
    cuffl_dir = sprintf('%scufflinks/drosophila/D_melanogaster_L3.%s/sanitized_out/', ...
                        data_dir, exp);
   case 'elegans',
    cuffl_dir = sprintf('%scufflinks/elegans/elegans.%s/sanitized_out/', ...
                        data_dir, exp);
   case 'thaliana',
    cuffl_dir = sprintf('%scufflinks/', data_dir);

  end
  cuffl_fn = 'transcripts.gtf';
  % load predictions
  cuffl = load([cuffl_dir cuffl_fn '.mat'], 'genes');
  % convert gene structure from closed intervals to half-open intervals
  cuffl.genes = closed_to_half_open(cuffl.genes);
  % load mTIM's config file
  load([pred_dir fname '.mat'], 'CFG');

if 0,
  % TODO
  keyboard
  train_time = '2011-04-13_13h10'; pred_time = '2011-04-28_14h56';
  exp_dir = 'mtim_exp_pair_04Apr';
  pred_dir = sprintf('%s%s.melanogaster.%s/galaxy/predictions/', data_dir, exp_dir, exp);
  fname = sprintf('genes_pred%s_trained%s', pred_time, train_time);
  cuffl = load([pred_dir fname '.mat'], 'genes');
  load([pred_dir fname '.mat'], 'CFG');
end

% hack for new version
CFG.xval_data_dir = CFG.out_dir;

  data_fn = CFG.xval_data_dir;
  idx = find(CFG.read_map_file == '/', 1, 'last');
  data_fn = [data_fn '/xval_fold1/test_data.mat'];
  load(data_fn, '-mat', 'test_chunks');
  fprintf('Using test chunks from %s\n', data_fn);

  if isfield('expr', pred.genes),
    pred.genes = norm_expr_levels(pred.genes,10);
  else
    for g=1:length(pred.genes),
      pred.genes(g).epxr = 10;
    end
  end

  if isfield('expr', cuffl.genes),
    cuffl.genes = norm_expr_levels(cuffl.genes,10);
  else
    for g=1:length(cuffl.genes),
      cuffl.genes(g).epxr = 10;
    end
  end

  % merge adjacent test chunks to save time during feature gathering
  MAX_LEN = 20000;
  c = size(test_chunks,1);
  while c > 1,
    if test_chunks(c,3)-test_chunks(c,2)+1 > MAX_LEN,
      c = c-1;
      continue
    end
  
    % merge if test chunks are on the same chromosome and adjacent
    if test_chunks(c,1) == test_chunks(c-1,1) ...
          && test_chunks(c,2)-1 == test_chunks(c-1,3),
      
      test_chunks(c-1,3) = test_chunks(c,3);
      % clear test chunk id
      test_chunks(c-1,4) = nan;
      test_chunks(c,:) = [];
    end
    c = c-1;
  end
  r = randperm(size(test_chunks,1));
  test_chunks = test_chunks(r,:);

  addpath ../data_preparation
  ag_int = [anno.genes.start; anno.genes.stop];
  pg_int = [pred.genes.start; pred.genes.stop];
  cg_int = [cuffl.genes.start; cuffl.genes.stop];
  for c=1:size(test_chunks,1),
    chunk = test_chunks(c,:);
    signal = fill_chunks(chunk, CFG);
    offset = chunk(2);
    ag = anno.genes([anno.genes.chr_num]==chunk(1) ...
                    & ag_int(2,:)>=chunk(2) & ag_int(1,:)<=chunk(3));
    pg = pred.genes([pred.genes.chr_num]==chunk(1) ...
                    & pg_int(2,:)>=chunk(2) & pg_int(1,:)<=chunk(3));
    cg = cuffl.genes([cuffl.genes.chr_num]==chunk(1) ...
                    & cg_int(2,:)>=chunk(2) & cg_int(1,:)<=chunk(3));    
    plot_region(signal, offset, ag, pg, cg)
    keyboard
  end

  keyboard
end



%%% Expression-level dependent evaluation
if LEVEL_EVAL,
  lvl_eval = struct();
  if isfield(pred.genes, 'expr'),
    pred.genes = norm_expr_levels(pred.genes, 10);
    levels = unique([pred.genes.expr]);
    
    %{
    for i=1:length(levels),
      l = levels(i);
      tmp_genes = pred.genes(round(gene_levels) == l);
      lvl_eval = struct();
      lvl_eval.intron = intron_eval(anno.genes, tmp_genes);
      fprintf('Intron evaluation for expression level %i:\n', l);
      lvl_eval.intron
    end
    %}

    for i=0:length(levels)-1,
      if i==0,
        l = 0;
      else
        l = levels(i);
      end
      tmp_genes = pred.genes(gene_levels > l);
      lvl_eval(i+1).intron = intron_eval(anno.genes, tmp_genes);

      lvl_eval(i+1).transcr = transcript_eval(anno.genes, tmp_genes);
      fprintf('Transcript evaluation for expression level %i-%i:\n', l+1,length(levels));
      lvl_eval(i+1).transcr

      lvl_eval(i+1).gene = gene_eval(anno.genes, tmp_genes);
    end
  end

  switch method,
   case 'mtim',
    glyph = 'p';
   case 'cufflinks',
    glyph = 'x';
   case 'scripture',
    glyph = '+';
  end
  
  col = [0.0 0.7 0.1; ...
         0.8 0.0 0.0;
         0.9 0.6 0.0];
  
  figure
  hold on
  
  in = [lvl_eval.intron];
  rec = [in.sens];
  prec = [in.prec];
  plot (rec, prec, glyph, 'Color', col(1,:), 'MarkerSize', 7)
  
  tr = [lvl_eval.transcr];
  rec = [tr.sens];
  prec = [tr.prec];
  plot (rec, prec, glyph, 'Color', col(2,:), 'MarkerSize', 7)
  
  
  gn = [lvl_eval.gene];
  rec = [gn.sens];
  prec = [gn.prec];
  plot (rec, prec, glyph, 'Color', col(3,:), 'MarkerSize', 7)
  
  title('Expression-dependant evaluation');
  xlabel('recall');
  ylabel('precision');
  grid on;
  axis([0 1 0 1]);
  
  fn_fig = sprintf('prec_rec_level_%s_%s.%s.eps', method, organism, exp);
  print('-depsc', [eval_dir fn_fig]);
  
  keyboard
end
