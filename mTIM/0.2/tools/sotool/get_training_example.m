function [Xi, Yi] = get_training_example(i, trainset)

idx = find(trainset.exm_id_intervals(:,1)==trainset.train_exm_ids(i));
idx = trainset.exm_id_intervals(idx,2):trainset.exm_id_intervals(idx,3);

Xi.obs_seq = trainset.signal(:,idx);
Yi.true_label_seq = trainset.label(idx);
Yi.true_state_seq = trainset.state_label(idx);
