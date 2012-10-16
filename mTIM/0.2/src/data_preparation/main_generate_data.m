function main_generate_data(CFG)
% Main script for data generation

% prepare label sequences
prepare_label(CFG);
fprintf('\n***************Label Prepared***************.\n');

% subdivide data into chunks
chunks = make_chunks(CFG);

fprintf('\n***************Made chunks***************.\n');
% assign chunks to cross-validation folds
split_data(CFG, chunks);

fprintf('\n***************Splitted data***************.\n');
% throw away some of the training chunks 
% BEFORE they are filled with label and feature information
filter_train_chunks(CFG);

fprintf('\n**************Filter trained***************.\n');
% fill chunks with label and feature sequences
generate_feature_data(CFG);

fprintf('\n**************Generated feature data***************.\n');

% eof
