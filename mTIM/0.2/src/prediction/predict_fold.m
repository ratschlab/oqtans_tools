function genes = predict_fold(predictor_dir, CFG)
% Predicts label sequences given a trained HM-SVM model
%
% predictor_dir -- a  directory where a trained HM-SVM model has been
%   saved
% returns a genes struct containing the predicted gene entries
%
% written by Georg Zeller & Gunnar Raetsch, MPI Tuebingen, Germany, 2008-2009
% written by Nico Goernitz, Berlin Institute of Technology, 2012
%   - migration from rproc to matgrid toolbox

% seed for random number generation
rand('seed', 11081979);

tic

LABELS = get_label_set_mTIM();

if predictor_dir(end) ~= '/',
  predictor_dir(end+1) = '/';
end

% load the final predictor
if isempty(CFG.iteration)
    % check if unconverged version exists
    fn = sprintf('%ssosvm_final.mat', predictor_dir);
else
  fn = sprintf('%ssosvm_iter%i.mat', predictor_dir, CFG.iteration);
end

% TODO: proper file check!
if ~exist(fn,'file'),
    load(fn(1:end-4), '-mat', 'PAR', 'information'); %'score_plifs', 'transition_scores');
else
    load(fn, '-mat', 'PAR', 'information');%'score_plifs', 'transition_scores');
end

model = information{end}.model;
score_plifs = model.score_plifs;
transition_scores = model.transition_scores;
fprintf('Loaded predictor from %s...\n', fn);

% TODO make this configurable
warning('symm_plifs is OFF.. TODO.');
%[score_plifs transition_scores] = symm_plifs(score_plifs, transition_scores, PAR);

if CFG.verbose>=2,
  disp(PAR);
  if CFG.verbose>=4 && ~CFG.grid_use,
    state_model = make_model_mTIM(PAR);
    view_model(state_model, score_plifs, PAR, transition_scores);
    keyboard
  end
end

if CFG.predict_genome,
  fprintf('Using equally-spaced genomic test chunks\n');
  chunks = make_genome_chunks(CFG);
else
  data_fn = PAR.data_file;
  idx = find(data_fn == '/', 1, 'last');

  if (CFG.predict_vald),
    data_fn = [data_fn(1:idx) 'vald_data.mat'];
    load(data_fn, '-mat', 'vald_chunks');
    fprintf('Using validation chunks from %s\n', data_fn);
    chunks = vald_chunks;
  else
    data_fn = [data_fn(1:idx) 'test_data.mat'];
    load(data_fn, '-mat', 'test_chunks');
    fprintf('Using test chunks from %s\n', data_fn);
    chunks = test_chunks;
  end

  warning('HACK: filter test chunks that have chunks(:,1) > CFG.num_chr.');
  inds = find(chunks(:,1) <= CFG.num_chr);
  chunks = chunks(inds,:);

  num_orig_chunks = size(chunks,1);
  % merge adjacent test chunks to save time during feature gathering
  MAX_LEN = 500000;
  c = size(chunks,1);
  while c > 1,
    if chunks(c,3)-chunks(c,2)+1 > MAX_LEN,
      c = c-1;
      continue
    end
  
    % merge if test chunks are on the same chromosome and adjacent
    if chunks(c,1) == chunks(c-1,1) ...
          && chunks(c,2)-1 == chunks(c-1,3),
      
      chunks(c-1,3) = chunks(c,3);
      % clear test chunk id
      chunks(c-1,4) = nan;
      chunks(c,:) = [];
    end
    c = c-1;
  end
end

if isfield(CFG, 'specific_strand'),
  load(fn, '-mat', 'state_model');
  ige_state = state_model(1);
  assert(isequal(ige_state.name, 'IGE'));
  assert(all(ige_state.trans_scores > 0));
  ige_succ = [state_model(ige_state.successors).label];
  switch CFG.specific_strand,
   case '+',
    idx = find(ige_succ == LABELS.exon_C);
    assert(length(idx) == PAR.num_levels);
    transition_scores(ige_state.trans_scores(idx)) = -inf;
   case '-',
    idx = find(ige_succ == LABELS.exon_W);
    assert(length(idx) == PAR.num_levels);
    transition_scores(ige_state.trans_scores(idx)) = -inf;    
   otherwise
    error('unknown strand %s', CFG.specific_strand)
  end
  fprintf('Making predictions for %s strand...\n\n', CFG.specific_strand);
else
  fprintf('Making predictions for both strands...\n\n');
end


% make predictions
fprintf('Predicting on data from %s\n', CFG.read_map_file);
fprintf('Processing %i chunks...\n', size(chunks,1));
for c=1:size(chunks,1),
  ARGS.chunk             = chunks(c,:);
  ARGS.CFG               = CFG;
  ARGS.transition_scores = transition_scores;
  ARGS.score_plifs       = score_plifs;
  ARGS.state_model       = model;
  ARGS.pred_PAR          = PAR;
 
  if CFG.grid_use,
      jobinfo(c) = mgsub({'pred_genes'}, 'predict_chunk', ...
        {ARGS}, 'qsub_opts', sprintf('-l h_vmem=%iG',CFG.grid_memreq));
      fprintf('  submitted job for chunk %i (%2.1f%%)\r',c, 100*c/size(chunks,1));
  else
    pred_genes{c} = predict_chunk(ARGS);
    if CFG.verbose >= 2,
      fprintf('  processed chunk %i (%2.1f%%)\r', c, 100*c/size(chunks,1));
    end
  end
end

if CFG.verbose >= 1,
  fprintf('  successfully processed %i chunks                      \n\n', c);
end

if CFG.grid_use,
  fprintf('Waiting for jobs to finish...');
  while true,
      if mgshow(jobinfo,1,15), break; end;
  end
  fprintf('Done!\n');
    
  fprintf('  collecting results\n');
  pred_genes = cell(1,length(jobinfo));
  for c=1:size(chunks,1),
    pred_genes{c} = mgresult(jobinfo(c));
  end
end

genes = [pred_genes{:}];
fprintf('  collected %i genes\n', length(genes));
fprintf('Whole-genome predictions took %.1f sec in total.\n\n', toc);
for g=1:length(genes);
  genes(g).name = sprintf('mTIM:gene_%i', g);
  genes(g).transcripts = {sprintf('mTIM:transcript_%i', g)};
end


% eof
