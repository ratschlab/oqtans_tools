function main_generate_data(CFG)
% Main script for data generation

% prepare label sequences
prepare_label(CFG);

% subdivide data into chunks
chunks = make_chunks(CFG);
% assign chunks to cross-validation folds
split_data(CFG, chunks);

% throw away some of the training chunks 
% BEFORE they are filled with label and feature information
filter_train_chunks(CFG);
% fill chunks with label and feature sequences
generate_feature_data(CFG);


% eof
