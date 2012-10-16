function prepare_data(dirname, organism)
% Writes example snippets into one big file according to the desired
% label format. Also collects some basic statistics.
%
% Expects the data examples in the directory 'dirname/organism/mat' 
% in the format 'example_%i.mat' where %i is a unique number.
%
%
%
% written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany, 2011

addpath('model');
LABELS = get_label_set();

% file prefix
file_prefix = 'example_';
fp_len = length(file_prefix);

% input directory where all examples lie
%base_dir = ['../../data/Escherichia_coli_BW2952_uid59391'];
%base_dir = ['../../data/Enterobacter_638_uid58727'];
base_dir = [dirname '/' organism];
data_dir = [base_dir '/mat/']
% output file name
data_output = [base_dir '/data.mat'];

% map m
% 0,1,2 start codons
% 61 62 63 stop codons
% m{codon_idx+1} = codon
fprintf('Load mapping file (expect idx_to_codon.mat in the current directory)..');
load(['idx_to_codon.mat']);
fprintf('Done!\n');

% get directory statistics
D = dir(data_dir);
fprintf('Data directory "%s" contains %i entries.\n',data_dir,length(D));

% output variables
exm_id = [];
label = [];
pos_id = [];
subset_id  =[];
signal = [];

% counters & statistics
max_exm_id = -inf;
min_exm_id = +inf;
num_exm = 0;

cnt_skip_exm = 0;

% integer values for start and stop codons
sCodonInds = [0 1 2 3 4];
tCodonInds = [61 62 63];

% load data files
for i=1:length(D),
  
  % load only valid files
  fname = D(i).name;
  
  if (length(fname)>=fp_len && strcmpi(fname(1:fp_len),file_prefix))
    % current example number is
    num = str2num(fname(fp_len+1:end-4));
  
    load([data_dir D(i).name]);
    interval = uint32(interval);
    codons_idx2 = uint8(codons_idx);
    codons_idx = zeros(1,length(sequence));
    codons_idx(1:length(codons_idx2)) = codons_idx2;
    % sequence: [1x1902 char] - die Primärsequenz
    % interval: [2x1 int64] - das Interval des CDS [start, stop[
    % mask: [1902x1 double] - das Label nochmal als 0, 1 Maske (1 --> CDS)
    % codons: [1900x3 char] - die codons an positionen i=0,...,n-2
    % codons_idx: [1900x1 int64] - die codons an positionen i=0,...,n-2, kodiert als integer
    
    % setup feature signal
    len = length(sequence);
    sig = zeros(1,len);
    sig(1:length(codons_idx)) = codons_idx';

    % CONVERT LABEL SEQUENCE
    % see model/get_label_set
    % LABELS.intergenic = -1;
    % LABELS.startCodon = 1;
    % LABELS.exonic = 2;
    % LABELS.stopCodon = 3;
    mask(mask==0) = LABELS.intergenic; % mark intergenic
    mask(mask==1) = LABELS.exonic; % mark exonic
    
    % find start and stop codon (strong assumption
    % that startCodon=ATG and stopCodon=AAT|GAT
    sind = interval(1)+1;
    tind = interval(2)+1;
    
    if (mod((tind-sind),3)~=0), 
      % skip example if not inframe
      cnt_skip_exm=cnt_skip_exm+1;
      continue; 
    end;

    if any(ismember(codons_idx(sind:3:tind-3), tCodonInds)), 
      % skip example if not inframe
      cnt_skip_exm=cnt_skip_exm+1;
      continue; 
    end;
    
    % skip example if anything is unexpected
    if (~any(codons_idx(sind)==sCodonInds)), cnt_skip_exm=cnt_skip_exm+1; continue; end;
    if (mask(sind)~=LABELS.exonic || mask(sind-1)~=LABELS.intergenic), cnt_skip_exm=cnt_skip_exm+1; continue; end;
    if (~any(codons_idx(tind)==tCodonInds)), cnt_skip_exm=cnt_skip_exm+1; continue; end;
    if (mask(tind)~=LABELS.intergenic || mask(tind-1)~=LABELS.exonic), cnt_skip_exm=cnt_skip_exm+1; continue; end;

    % set start and stop labels
    mask(sind) = LABELS.startCodon;
    mask(tind) = LABELS.stopCodon;
    
    % set id
    id = num*ones(1,len);
    
    % add to the pool
    signal = [signal sig];
    exm_id = [exm_id id];
    label = [label mask'];
    
    % increase counters
    num_exm = num_exm+1;
    if (num>max_exm_id), max_exm_id=num; end;
    if (num<min_exm_id), min_exm_id=num; end;
  end
end

% build example intervals
all_exm_ids = unique(exm_id);
fprintf('Building intervals array..');
exm_id_intervals = zeros(length(all_exm_ids),3);
for i=1:length(all_exm_ids),
  idx = find(exm_id==all_exm_ids(i));
  exm_id_intervals(i,:) = [all_exm_ids(i), idx(1), idx(end)];
end
fprintf('Done!\n');

% plot some statistics
fprintf('%i examples converted.\n',num_exm);
fprintf('%i examples skipped.\n',cnt_skip_exm);
fprintf('Min/Max exm id was %i/%i.\n',min_exm_id,max_exm_id);
fprintf('Length of the sequence is %i.\n',length(label));


% empty signals
pos_id = zeros(1,length(label));
subset_id = ones(1,length(label));

% save data
fprintf('Write output to "%s".\n',data_output);
save(data_output,'exm_id','label','signal','pos_id','subset_id','exm_id_intervals');

% plot some distribution statistics
start_feats = logical(zeros(1, length(label)));
stop_feats  = logical(zeros(1, length(label)));
codon_feats = logical(zeros(1, length(label)));
% out-of-frame
oof_feats   = logical(zeros(1, length(label)));

s_idx = find(label==LABELS.startCodon);
t_idx = find(label==LABELS.stopCodon);
for f=1:length(s_idx),
  start_feats(s_idx(f)) = 1;
  for p=s_idx(f)+3:3:length(label),
    if label(p)==LABELS.stopCodon,
      stop_feats(p) = 1;
      break
    end
    assert(p < t_idx(f));
    codon_feats(p) = 1;
    oof_feats(p+1) = 1;
    oof_feats(p+2) = 1;
  end
end

assert(~any(start_feats & stop_feats & codon_feats & oof_feats));
% intergenic
ige_feats = logical(ones(1, length(label)));
ige_feats(start_feats | stop_feats | codon_feats | oof_feats) = 0;


start_feats = signal(start_feats);
stop_feats  = signal(stop_feats);
codon_feats = signal(codon_feats);
oof_feats   = signal(oof_feats);
ige_feats   = signal(ige_feats);

n_start = hist(start_feats, 0:63);
n_stop  = hist(stop_feats, 0:63);
n_codon = hist(codon_feats, 0:63);
n_oof   = hist(oof_feats, 0:63);
n_ige   = hist(ige_feats, 0:63);

figure(2)
subplot(5,1,1);
bar(0:63, n_start);
title('Start Codon State')
subplot(5,1,2);
bar(0:63, n_stop);
title('Stop Codon State')
subplot(5,1,3);
bar(0:63, n_codon);
title('In-Frame Codon States')
subplot(5,1,4);
bar(0:63, n_oof);
title('Out-of-Frame Exon States')
subplot(5,1,5);
bar(0:63, n_ige);
title('Intergenic State')


