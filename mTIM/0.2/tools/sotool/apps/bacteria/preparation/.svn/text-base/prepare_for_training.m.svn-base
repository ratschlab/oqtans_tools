function prepare_for_training(reldir, num_reps, num_vald, num_test)
% Loads a 'cds_genes.mat' file in directory 'reldir' and
% prepares the data of the organisms for training of the sosvm.
%
% The first step ...
%
% The second step prepares the training/validation and testing
% slices. Here, it is crucial to take COGs into account:
% there shouldn't be genes within one slice (also run) which
% have identical COGs either within one organism or across.
%
% Outputs:
%   runs                : [num_reps x num_orgs x (num_train+num_vald+num_test)]
%   exm_id_intervals    : [2 x #total_examples]
%   signal              : [#dim x all_sequences] 
%   label               : [1 x all_sequences]
%   state_label         : [1 x all_sequences]
%
%
% author nico goernitz, MPI Tuebingen, TU Berlin, 2011

addpath('../model');

% 1. CONSTANTS
fname = 'cds_genes.mat';
reldir = '../../../out/tax4/';
load([reldir fname]);

num_reps = 10;
num_train = 40;
num_vald = 100;
num_test = 200;


% 2. CONVERSION/EXTENSION OF LOADED VARIABLES
% extend the taxonomy with single train on leafs
num_orgs = size(org_interval,1);
tax_idx_names = regexp(org_names, '([\t]+)', 'split');
tax_idx_names = tax_idx_names(1:end-1);
org_idx_names = tax_idx_names((end-num_orgs+1):end);

% data type conversions
sequence = int8(sequence);
interval = uint32(interval);
relative_interval = uint32(relative_interval);

% fix example_id bug
example_id = 1:length(example_id);

% fix interval & org_interval bug
add = 0;
last_end_pos = 0;
for i=1:length(example_id),
    if (last_end_pos>(interval(i,1)+add)),
        % fix broken interval 
        add = last_end_pos - interval(i,1) + 1;
        fprintf('Fix broken interval last_pos=%i new_pos=%i.\n',last_end_pos,interval(i,1)+add);
    end
    interval(i,1) = interval(i,1) + add;
    interval(i,2) = interval(i,2) + add;
    last_end_pos = interval(i,2);
end
assert(interval(end,2)==length(sequence))



% 3. CONVERT/CHECK SEQUENCE
fprintf('\nCheck and convert sequence.\n');

assert(length(start_codon) == length(stop_codon));
assert(length(start_codon) == length(example_id));


% 3.1. filter all examples that are somehow corrupted
% and cut sequences (-2 because handling without codon information is not implemented yet)

exclude_ids = [];
% for all genes 
for i=1:length(example_id),
    start_pos = interval(i,1)+relative_interval(i,1);
    stop_pos = interval(i,1)+relative_interval(i,2);

    % start and stop should be some nukleotides apart
    if (stop_pos < (start_pos+4)),
        % remove
        exclude_ids = [exclude_ids, example_id(i)];
    end
    
    % exon should not start at the first two positions
    if (start_pos <= interval(i,1)+3),
        % remove
        exclude_ids = [exclude_ids, example_id(i)];
    end

    % exon has to end some positions before end of sequence -2
    if (stop_pos >= interval(i,2)-5),
        % remove
        exclude_ids = [exclude_ids, example_id(i)];
    end

    % check if inframe length is mod 3
    if (mod((stop_pos-start_pos),3)~=0),
        % remove
        exclude_ids = [exclude_ids, example_id(i)];
    end
end
exclude_ids = unique(exclude_ids);
fprintf('%i non-usable genes found (which are excluded from training/vald and testing).\n',length(exclude_ids));

% 3.2. create label sequence
% setup signal(observation vector) -2 positions not covered by 3-mers
signal = zeros(1,length(sequence) - length(example_id)*2);
signal = uint8(signal);

% example boundaries
exm_id_intervals = zeros(2,length(example_id));

LABELS = get_label_set();
label = LABELS.intergenic*ones(1,length(signal));
label = int8(label);

% for all genes 
last_end_pos = 0;
for i=1:length(example_id),

    % setup signal vector which is mainly the codons_idx vector
    % but cut last 2 positions of each example sequence
    len = interval(i,2)-interval(i,1)-2;
    exm_id_intervals(1,i) = last_end_pos + 1;
    exm_id_intervals(2,i) = exm_id_intervals(1,i) + len;
 
    signal(exm_id_intervals(1,i):exm_id_intervals(2,i)) = codons_idx(interval(i,1):interval(i,2)-2);
    last_end_pos = exm_id_intervals(2,i);

    % only label useful sequences
    if (sum(exclude_ids==example_id(i))==0),
        start_pos = exm_id_intervals(1,i)+relative_interval(i,1);
        stop_pos = exm_id_intervals(1,i)+relative_interval(i,2);

        label(start_pos:stop_pos) = LABELS.exonic;
        label(start_pos) = LABELS.startCodon;
        label(stop_pos) = LABELS.stopCodon;

        % check if the start and stop codons are valid
        assert( any(start_codons_idx==signal(start_pos)) );
        assert( any(stop_codons_idx==signal(stop_pos)) );

        assert( mod(stop_pos-start_pos,3)==0 );
    end
end
assert(sum(label==LABELS.startCodon)==sum(label==LABELS.stopCodon));
assert(length(signal)==exm_id_intervals(2,end));

% 3.3. create state sequence
fprintf('Converting label to state sequences..');
PAR.num_features = 1;
[state_label, problems] = labels_to_states(label, make_model(PAR), signal, PAR);
assert(isempty(problems));
fprintf('Done!\n');




% 4. PREPARE TRAIN/VALD/TEST SETTING
fprintf('\nCreate train/vald/test runs.\n');
samples = num_train+num_vald+num_test;

% the train/vald/test ids for all runs
runs = -ones(num_reps, num_orgs, samples);

% a corresponding tensor for homology checking
cogs = {};

org_exms = zeros(1,example_id);
for o=0:(num_orgs-1),
    org_exms(o+1) = sum(org_id==o);
end
[foo, org_order] = sort(org_exms);
org_order = org_order-1;

for r = 1:num_reps,
    cogs{r} = [];
    fprintf('Draw training/validation and testing set for run %i.\n',r);
    % heuristic: sort the organisms with respect to their amount of
    % examples    
    for o = org_order,
        % org id starts with 0

        % get all the example ids belonging to the current organism
        exm_ids = example_id(org_id==o);
        cog_ids = cog(org_id==o);
        total_exms = length(exm_ids);
    
        % take care of excluded examples
        bitmask = ~ismember(exm_ids,exclude_ids);
        exm_ids = exm_ids(bitmask==1);
        cog_ids = cog_ids(bitmask==1);
        incl_exms = length(exm_ids);

        % remove all examples which have a cog id that is already used
        bitmask = ~ismember(cog_ids, cogs{r});
        bitmask = bitmask | [cog_ids==-1];
        exm_ids = exm_ids(bitmask==1);
        cog_ids = cog_ids(bitmask==1);
        filter1_exms = length(exm_ids);

        % remove all duplicates within cog_ids
        [foo, inds] = unique(cog_ids);
        inds = [inds; find(cog_ids==-1)];
        % remove double indices
        inds = unique(inds);
        exm_ids = exm_ids(inds);
        cog_ids = cog_ids(inds);
        filter2_exms = length(exm_ids);

        % test if there are still enough examples
        if (length(exm_ids) < samples),
            fprintf('Not enough training examples (%i/%i).\n',length(exm_ids),samples);
            fprintf('Restarting the skript may help.\n');
        end
        assert(length(exm_ids) >= samples);

        % shuffle and cut sample set data
        rndinds = randperm(length(exm_ids));
        exm_ids = exm_ids(rndinds(1:samples));
        cog_ids = cog_ids(rndinds(1:samples));

        fprintf('Org(%i): total/valid/filter1/filter2/samples = %i/%i/%i/%i/%i.\n',o,total_exms,incl_exms,filter1_exms,filter2_exms,length(exm_ids));

        runs(r,o+1,:) = exm_ids;
        cogs{r} = [cogs{r}, cog_ids];
    end
end



% 5. SAVE 
fname = sprintf('data_reps%i_tr%i_val%i_tst%i.mat',num_reps,num_train,num_vald,num_test);
save([reldir fname], 'taxonomy', 'tax_idx_names', 'org_idx_names', 'cogs','runs', 'label','state_label','signal',...
    'exm_id_intervals','num_train','num_vald','num_test','example_id','org_id','exclude_ids');


