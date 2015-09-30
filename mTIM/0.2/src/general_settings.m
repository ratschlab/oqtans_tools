function CFG = general_settings(...
            exp_name, time_stamp, ...
            create_prep_dirs, create_exp_dirs)
% In:
%   exp_name  : (String) name of the experiment (e.g. 'elegans.6.filtered');
%               !effects directory names and the 
%               specialized settings file: %exp_name%_settings.m
%   time_stamp: corresponds to a CFG.out_train_dir
%               if not set OR EMPTY the config will use the current date
%   create_prep_dirs: create data preparation dirs if set ~0
%   create_exp_dirs : create experiment dirs if set ~0
%
% Out:
%   CFG       : (Struct) a structure containing all necessary information
%               for any step of an experiment
%
%
% Sets ALL paths and general settings 
% in one location for all stages of an experiment:
%
%  1) data preparation
%  2) training/model selection
%  3) prediction
%  4) evaluation
% 
% 2) and 3) can be done either sequentially or parallel
% 

% locate src directory based on the location of this script
this_name = mfilename();
this_path = mfilename('fullpath');
idx = find(this_path=='/', 1, 'last');
src_dir = this_path(1:idx);
assert(exist(src_dir,'dir')>0);

% if the experiment name starts with an
% absolut path, then use the this path
% as main directory and the last part (sub-sequence
% after the last '/') is the exp_name
%if (exp_name(1)=='/'),
%    idx = find(exp_name=='/', 1, 'last');
%    main_dir = exp_name(1:idx);
%   exp_name = exp_name(idx+1:end);
%    fprintf('Absolut path detected:\n  Using %s as main_dir and\n  %s as exp_name.\n',main_dir,exp_name);
%else
%    idx = find(this_path=='/', 2, 'last');
%    main_dir = this_path(1:idx);
%    fprintf('Relative path detected:\n  Using %s as main_dir.\n',main_dir);
%end

%out_dir = [main_dir 'out/'];
%if (~exist(out_dir,'dir')), mkdir(out_dir); end;
%out_dir = [main_dir 'out/' exp_name '/'];
%if (~exist(out_dir,'dir')), mkdir(out_dir); end;

CFG.out_dir = exp_name;
CFG.exp_name = exp_name;

if (~exist('time_stamp','var') || isempty(time_stamp)), 
    time_stamp = datestr(now,'yymmdd_HHMM'); 
end
CFG.start_time = time_stamp;

CFG.out_train_dir = [exp_name '/' CFG.start_time '/'];
CFG.out_pred_dir = CFG.out_train_dir;
CFG.pred_filename = 'prediction.mat';
CFG.vald_filename = 'validation.mat';

%fprintf('Main directory is %s.\n',main_dir);
fprintf('Out directory is %s.\n',exp_name);
fprintf('Source directory is %s.\n',src_dir);


% add subdirs
model_dir = [src_dir 'model'];
common_dir = [src_dir 'utils'];
pred_dir = [src_dir 'prediction'];
prep_dir = [src_dir 'data_preparation'];
hmsvm_dir = [src_dir 'sotool'];

addpath(model_dir);
addpath(common_dir);
addpath([src_dir 'utils2']);
addpath(prep_dir);
addpath([src_dir 'training']);
addpath(pred_dir);
addpath([src_dir 'evaluation']);
%addpath([src_dir 'rproc']);
addpath(hmsvm_dir);

% paths needed during PARALLEL training/prediction
CFG.PAR.include_paths = {src_dir, model_dir, prep_dir, common_dir, hmsvm_dir};

% ugly hack!: if the current training directory already exists 
% then this is a prediction task (only needed for rproc_start_dir)
is_pred_task = 0;
if (~create_prep_dirs && ~create_exp_dirs), is_pred_task = 1; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRAINING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% cross-validation / model selection parameters
% names of parameters to be independently specified (possibly differently
% between training runs) 
CFG.param_names = {'C_small', ...
                   'C_smooth', ...
                   'C_coupling', ...
                   'num_train_exm', ...
                   'reg_type', ...
                   'level_loss' ...
                  };
% parameter combinations to be used for independent training
% (typically overwritten)
CFG.train_params = { ...
  [0.01],   [0.5],   [0.05],  [20], 'QP', [1]; ... 
};

CFG.num_vald_exm = 0;

% setup model directories 
CFG.model_dirs = {};

rand('seed',11081979);

% train the following classifier
CFG.PAR.classifier = @train_primal_sosvm;

% bundle method settings
CFG.PAR.bmrm.useMonotConstr = 0;  % use monotonicity constraints (0=no 1=yes)
CFG.PAR.bmrm.eps = 1e-1;          % relative gap
CFG.PAR.bmrm.useAggregation = 1;  % use aggregation (set appropriate K)
CFG.PAR.bmrm.K = 168;             % number of cutting planes used -1 + aggregated plane
CFG.PAR.bmrm.checkRealObjMod = 5; % heuristic 1
CFG.PAR.bmrm.monotMargin = 0.01;  % heuristic 2
% bmrm linesearch settings
CFG.PAR.bmrm.useLs = 1;           % use linesearch (0=no 1=yes)
CFG.PAR.bmrm.lsEps = 1e-3;        % break if parabola approx gap is smaller than lsEps
CFG.PAR.bmrm.lsSteps = 50;        % max. number of iterations
CFG.PAR.bmrm.lsTheta = 0.25;      % w_new = (1-lsTheta)*w_ls + lsTheta*w_old;


% number of discrete expression levels
CFG.PAR.num_levels = 5; % something like 3-10;
% number of supporting points for each scoring function
CFG.PAR.num_plif_nodes = 2*CFG.PAR.num_levels; %CFG.PAR.num_levels+2 %2*CFG.PAR.num_levels;

% include mate-pair information for paired reads as a learning features
CFG.PAR.use_pair_feats = 0;
% include data on genomic repeats as a learning features
CFG.PAR.use_repeat_feats = 0;
% use splice feats
CFG.PAR.use_splice_feats = 1;
% use original filtered intron span feature
CFG.PAR.use_filtered_intron_feats = 1;
% use the binned intron span features
CFG.PAR.use_binned_span_feats = 0;
% use the cufflinks feature
CFG.PAR.use_cuffl_feats = 0;

% enforce the monotonicity for some feature scoring functions
% note that this can make the training problem more difficult
% (i.e. time-consuming) to solve!
CFG.PAR.enf_monot_score_funcs = 0;
% enforce monotonicity only in IGE state
CFG.PAR.enf_monot_score_funcs_IGE = 0;

% remove all positions where low_cover_blocks feature < PAR.gene_states_low_cover_cutoff
% before training (see train_hmsvm.m)
PAR.remove_low_cover_blocks = 0;

% use heuristic training procedure
CFG.PAR.constraint_margin = 5;
% options for enforcing a certain  state given a certain range of feature values
CFG.PAR.switch_features = 1;
if CFG.PAR.switch_features,
  % the maximum allowed value for the low-coverage block feature in exon
  % and intron states (forcing to decode these as intergenic regions)
  CFG.PAR.gene_states_low_cover_cutoff = 10; % has to be CFG.PAR because it 
                                             % is used in train_hmsvm as well
                                         % a value of 10 corresponds to a block
                                         % of ~2kb with mean coverage of 1 or 
                                         % a 1kb-region with 0-coverage
  % restrict the range of splice site predictions to be > -inf in the IF and
  % IL states to effectively enforce the canonical splice site dinucleotide
  % consensus
  CFG.PAR.enforce_splice_site_consensus = 1;
  % compute two matrices states x features specifying a range of allowed
  % feature values; if a feature value outside this range is encountered,
  % decoding the corresponding state will not be possible. Hence, too many
  % / too strict ranges can make decoding (& training) impossible!
  CFG.PAR.perm_feature_ranges = set_permitted_feature_ranges(CFG);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NUMBER of subsets used for cross-validation
CFG.num_xval_folds = 3;

% build xval dir names (for training and data-preparation!)
CFG.xval_dirs = {};
CFG.xval_train_dirs = {};

CFG.num_data_subsets = 3;

CFG.max_train_chunk_len = 250000;
CFG.fill_vald_chunks = 1;
CFG.fill_test_chunks = 1;

% genome-wide positional labelings for mTIM
CFG.label_fn = sprintf('%s/label', exp_name);

% parameters for generating the intron span feature
CFG.read_intron_span_max_intron_len = 200000; % TODO sufficient for worm&fly
CFG.read_intron_span_min_score = 25; % C. elegans & D. melanogaster
CFG.read_intron_span_max_mismatches = 0; 

% parameters for extracting splice site features from spliced read alignments
CFG.read_splice_site_max_intron_len = 200000; % TODO sufficient for worm&fly
CFG.read_splice_site_min_score = 10; 
CFG.read_splice_site_max_mismatches = 0; 

% parameters for retrieving read pair features from read alignments
CFG.read_pair_cover_max_intron_len = 200000;
CFG.read_pair_cover_min_score = 25;
CFG.read_pair_cover_max_mismatches = 0;

% threshold on total read coverage used to define low-coverage regions
CFG.low_cover_threshold = 5;  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREDICTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predict on test (rather than validation) set?
CFG.use_vald_data = 0;
% make test predictions for the whole genome?
CFG.predict_genome = 0;
% make predictions for a single strand only?
%CFG.specific_strand = '-';

% specific training time (corresponds to some CFG.out_train_dir)
CFG.train_time = time_stamp; 

% optimal model of mTIM to be selected for prediction (usually 1, unless
% several models were trained in the same run_crossvalidation.m script)
CFG.opt_model = 1;
% training iterations (take converged predictor if empty)
CFG.iteration = [];


%%% cluster submission vs. local computation
% option for using rproc tools to submit training jobs to a cluster
CFG.use_rproc = 0;
CFG.rproc_memreq         = 7000;
CFG.rproc_par.force_matlab = 1;
CFG.rproc_par.identifier = sprintf('mTIM_%s', CFG.exp_name);
% ugly hack 
CFG.rproc_par.start_dir  = [hmsvm_dir];
if (is_pred_task),  CFG.rproc_par.start_dir = [pred_dir]; end
CFG.rproc_par.verbosity  = 0;
CFG.rproc_par.priority   = 505;
CFG.rproc_time           = 24*60; % mins

%%% behavior during training
CFG.PAR.extra_checks = 0;
CFG.PAR.verbose      = 1;
CFG.PAR.check_acc    = 0;
CFG.PAR.submit_train = 1;
CFG.PAR.submit_vald  = 0;
CFG.PAR.submit_batch = 0;

CFG.verbose = 1;

% call the specialized settings file according to 'exp_name'
special_settings = [exp_name '_settings'];
if (exist(special_settings,'file')),
    fprintf('Calling special settings file: %s\n',special_settings);
    set_special_settings = str2func(special_settings);
    CFG = set_special_settings(CFG);
else
    warning('Could not find special settings file: %s\n',special_settings);
end

% create the directory structure for the current experiment
CFG = create_dirs(CFG,create_prep_dirs,create_exp_dirs);

% eof
